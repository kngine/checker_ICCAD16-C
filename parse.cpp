#include <cstdio>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#include "header.h"

extern int read_gds(const char *gds_file_path, layout *layoutLib);
//read gds file to layout

extern void clear_layers();
//clear the layer info of gds

extern bool isSimilar(clip *representative, clip *clip_x);
//return whether clip_x is similar to the representative

using namespace std;
using namespace boost;

#define MAXSIZE 10000 //max number of clusters
int ISVALID = true;

bool pointInPoly(polygon* p, dpoint pt)
{
    return pt.x >= p->ll.x && pt.y >= p->ll.y 
        && pt.x <= p->ur.x && pt.y <= p->ur.y;
}

bool cmpCluster(const vector<int> &a, const vector<int> &b)
{
    return a.size() > b.size() || (a.size() == b.size() && a[0] < b[0]);
}

void printClusters(vector<vector<int> > clusters)
{
    printf("=printing clusters:=\n");
    int cnt = 0;
    for(unsigned int i=0; i<clusters.size(); i++)
    {
        printf("cluster %d: ", i+1);
        cnt += clusters[i].size();
        sort(clusters[i].begin()+1, clusters[i].end());
        for(unsigned int j=0; j<clusters[i].size(); j++)
        {
            printf("%d ", clusters[i][j]);
        }
        printf("\n");
    }
    printf("Total clip #: %d\n", cnt);
}

void handleError(string s, int type=0)
{
    std::cout<<s;
    ISVALID = false;
}

void checkClusters(vector<vector<int> > &clusters, int repID[])
{
    printf("=checking cluster similarity...=\n");
    for(unsigned i=0; i<clusters.size(); i++)
    {
        printf("- now checking cluster %d  with representative clip %d -\n", i+1, repID[i]);
        clip *R = database.clips+repID[i];
        clip *A;
        for(unsigned int j=0; j<clusters[i].size(); j++)
        {
            A = database.clips+clusters[i][j];

            if(!isSimilar(R,A))//check whether similar according to the specified limit and in the specified mode (ACC/ECC)
            {
                char ss[256];
                sprintf(ss, "ERROR: similarity is not within the limit for clip %d\n", clusters[i][j]);
                handleError(ss); 
            }
        }
    }
    printf("=clusters checking finished=\n");
} 

int main(int argc, char **argv)
{
    /*------------------------------
    get arguments
    -------------------------------*/
    if(argc < 7)
    {
        printf("Usage: %s <input_gds_file> <ACC/ECC> <threshold> <clip_width> <clip_height> <overlay_gds_file> [representative_gds_file]\n", argv[0]);
        return -1;
    }
    database.mode=ACC;
    database.accThres=1;
    database.eccThres=0;

    char layout_file[100];
    char overlay_file[100];
    char rep_file[100]={0};
    strcpy(layout_file, argv[1]);
    
    if(strcmp(argv[2], "ACC")==0)
        database.mode = ACC;
    else if(strcmp(argv[2], "ECC")==0)
        database.mode = ECC;
    else
    {
        printf("\"%s\" in command line is not recognized\n", argv[2]);
        exit(-1);
    }

    if(database.mode==ACC)
        database.accThres = atof(argv[3]);
    else if(database.mode==ECC)
        database.eccThres = atof(argv[3]);

    database.clipW = atof(argv[4]);
    database.clipH = atof(argv[5]);

    strcpy(overlay_file, argv[6]);

    if(argc>=8) //representative_gds_file is optional
    {
        strcpy(rep_file, argv[7]);
    }

    /*------------------------------
    read file
    -------------------------------*/
    printf("parsing cluster file...\n");

    layout layoutLib[2]; //each layout in a layout lib is one layer
    //layoutLib[0] stores all the features, layoutLib[1] stores all the markers
    read_gds(layout_file, layoutLib);

    for (unsigned int i=0;i<layoutLib[1].polygons.size();i++)
    {
        layoutLib[1].polygons[i]->getCenter();
        layoutLib[1].polygons[i]->getLL();
        layoutLib[1].polygons[i]->getUR();
    }

    layoutLib[1].centerYX_sort();
    database.initClips(layoutLib[1].polygons.size());

    _extClip(layoutLib); //init data

    //a layout contains a set of polygons
    layout clusters[MAXSIZE]; //clusters[i] stores the layout of i-th cluster 
    layout rep[MAXSIZE]; //rep[i] stores the layout of i-th representative 
    
    clear_layers();
    read_gds(overlay_file, clusters);
    if(rep_file[0]!=0) 
    {
        read_gds(rep_file, rep); 
        //NOTE: this function will use layer ID to match rep_file and over_lay file
        //i.e., layer n of rep_file should contain the representative for layer n of overlay_file 
    }

    /*------------------------------
    process clusters and representatives
    -------------------------------*/
    const unsigned nMarker = layoutLib[1].polygons.size();
    unsigned nCLip = 0;

    typedef adjacency_list<vecS, vecS, undirectedS> my_graph; 
    my_graph g(nMarker * 2); //graph for matching
    std::vector<graph_traits<my_graph>::vertex_descriptor> mate(nMarker*2);  

    vector<vector<int> > clusterID; //clusterID[i] stores clip-IDs of the i-th cluster, clusterID.size() is #cluster
    int repID[MAXSIZE]; //repID[i] is the clip-ID of the representative of i-th cluster
    vector<polygon*> allClipPoly; //a vector store all the polygons of the clips in the overlay file
    
    for(unsigned int i=0; i<MAXSIZE; i++) //process the overlay file layer by layer, each layer is a cluster
    {
        if(clusters[i].polygons.size()==0) break; //all layers have been processed
        
        repID[i] = -1;
        if(rep_file[0]!=0) //if rep_file is set by command line
        {
            rep[i].polygons[0]->getLL();
            rep[i].polygons[0]->getUR();
        }
        vector<int> cluster0;
        for(unsigned int j=0; j<clusters[i].polygons.size(); j++) //process the polygons of i-th cluster
        {
            cluster0.push_back(nCLip); //record the clip-ID of current clip
            allClipPoly.push_back(clusters[i].polygons[j]);
            polygon *p1 = clusters[i].polygons[j]; //this is j-th polygon of i-th cluster
            p1->getLL();
            p1->getUR();
            p1->getCenter();
            if((p1->ur.x - p1->ll.x)!=database.clipW || (p1->ur.y - p1->ll.y)!=database.clipH)
            {
                handleError("ERROR: incorrect size for the following clip\n");
                p1->printPolygon();
            }

             if(repID[i]==-1)
            {
                if(rep_file[0]!=0)
                {
                    if(rep[i].polygons.size()>0)
                    {
                        polygon *p2 = rep[i].polygons[0];
                        if(p1->ll==p2->ll && p1->ur==p2->ur)
                        {
                            repID[i] = nCLip;
                        } 
                    }
                }
            }

            unsigned k;
            bool found = false;
            for (k=0;k<layoutLib[1].polygons.size();k++)
            {
                if(pointInPoly(layoutLib[1].polygons[k], clusters[i].polygons[j]->center))
                //clip's center is within the k-th marker    
                {
                    add_edge(k, nCLip + nMarker, g);
                    //g has 2*nClip vertices, the first half are marker-IDs, the second half are clip-IDs
                    //k-th marker may be the marker for nClip-th clip, so an edge is added to g
                    //NOTE: one clip's center may be within multiple markers, considering abutting/overlapping markers
                    //so a graph matching will be done to decide which is corresponding to which
                    //e.g.  g has 3 edges  (marker 0)--(clip-2) (marker 1)--(clip-2) (marker 1)--(clip-3)
                    //matching result will be (marker 0)--(clip-2) (marker 1)--(clip-3)

                    found = true;
                }
            }
            if(!found)
            {
                handleError("ERROR: the following clip's center is not within any marker\n");
                clusters[i].polygons[j]->printPolygon();
                // exit(-1);
            }
            nCLip++;
        }
        clusterID.push_back(cluster0); 
        if(repID[i]==-1) repID[i] = cluster0[0]; //if repID[i] is still not set, set it as the first clip encontered
    }

    if(nCLip!=nMarker) //these two numbers must be equal
    {
        handleError("ERROR: #clip!=#marker");
        exit(-1);
    }

    edmonds_maximum_cardinality_matching(g, &mate[0]);//ideally, the maximum matching is a perfect graph matching

    if(matching_size(g, &mate[0])!=nMarker)
    {
        handleError("ERROR: clips cannot be one-to-one maped to markers");
        exit(-1);
    }

    graph_traits<my_graph>::vertex_iterator vi, vi_end;
    for(tie(vi,vi_end) = vertices(g); *vi < nMarker; ++vi)
    {
        //(*vi)-th clip is mapped to (mate[*vi]-nMarker)-th clip
        database.clips[*vi].setBoder(allClipPoly[mate[*vi]-nMarker]);
        //extract the real layout (features) for the clip
        database.clips[*vi].cutClip();
        if(database.clips[*vi].polygons.size()==0)
        {
            printf("WARNING: clip %d is blank\n", (int)*vi);
        }
        // std::cout << "{" << *vi << ", " << mate[*vi]-nMarker << "}" << std::endl;
    }

    /*------------------------------
    print and check clusters
    -------------------------------*/
    sort(clusterID.begin(), clusterID.end(), cmpCluster);//sort by cluster size
    for(unsigned i=0; i<clusterID.size(); i++)
    {
        //according to matching result, map clip-IDs to marker-IDs
        for(unsigned j=0; j<clusterID[i].size(); j++)
        {
            clusterID[i][j] = mate[nMarker + clusterID[i][j]];
        }
        repID[i] = mate[nMarker + repID[i]];
        sort(clusterID[i].begin(), clusterID[i].end());
    }

    printClusters(clusterID);

    checkClusters(clusterID, repID);//check similarity for each cluster

    if(ISVALID)
        printf("Solution is valid.\n");
    else
        printf("!!!!!!Solution is invalid.!!!!!!\n");
    printf("----------Number of cluster: %d------------\n", (int)clusterID.size());

    return 0;
}
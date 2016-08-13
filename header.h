#include <cstdlib>
#include <vector>
#include <functional>
#include <algorithm>

class point
{
public:
    int x, y;
    unsigned char outdir;

    point();
    point(int x, int y, unsigned char dir = 0);
    bool operator  < (const point &T) const;
    bool operator == (const point &T) const;
    bool operator <= (const point &T) const;
    bool operator != (const point &T) const;
    point operator + (const point &T) const;
    point operator - (const point &T) const;
};

class dpoint
{
public:
    double x, y;
    bool operator  < (const dpoint &T) const;
    bool operator == (const dpoint &T) const;
    bool operator <= (const dpoint &T) const;
    bool operator != (const dpoint &T) const;
};

struct rect;

union addr_area
{
    rect* ref;
    int overlap_area;
};

struct rect
{
    point ll, ur;
	point reserve;
    long long waste;
    addr_area ext;
    std::vector<rect*> *upper, *lower;
    rect();
    rect(const point &ll, const point &ur);
};

class polygon
{
public:
    std::vector<point> pt_list;
    point ll, ur; // the lower-left/upper-right of the polygon
    dpoint center;
    int ID;       // the ID number of a polygon in the whole layout

    void getLL(); // get the lower left coordinate of the polygon. 
    void getUR(); // get the upper right coordinate of the polygon. 
    void getCenter(); // get center of the polygon. 
    void coordRevise();
    void getID(int n);
    void printPolygon()
    {
        printf("printing polygon :");
        int m = this->pt_list.size();
        for (int j=0;j<m;j++)
            printf("(%d,%d)", this->pt_list[j].x, this->pt_list[j].y);
        printf("\n");
    }
    ~polygon();
};

/*------------------------------
a layout is a set of polygons
-------------------------------*/
class layout
{
public:
    short layer;
    std::vector<polygon *>  polygons;
    point ll, ur;     // the lower-left and upper-right of the whole layout

    void initLayoutllur();
    void printPolygons();
    void llXY_sort(); // sort x first ,if x is the same then sort y.
    void llYX_sort();
    void centerYX_sort();
    ~layout();
};


/*--------------------------------------
a clip is a set of polygons and a marker
----------------------------------------*/
class clip
{
public:
    rect bBorder;                       // the rectangle information of the area of all clips within the marker  
    rect border;                        // the rectangle information of the clip
    rect marker;                        // the rectangle information of the marker
    point center;
    int markerType;                     // the type of the marker defined by NCS rule
    int numPoint;                       // the number of the sampling point

    std::vector<int> sampleID;
    std::vector<polygon *> polygons;    // store all the polygons intersected with the current small clip
    std::vector<polygon *> bPolygons;   // store all the polygons intersected with the big clip (only inside)
    std::vector<polygon *> oPolygons;   // store all the polygons intersected with the big clip (including outside)
    std::vector<rect *> rects;           
    
    void cutClip();
    void cutClip(layout*);
    void getBoder();
    void setBoder(polygon *);
    ~clip();
};

class Node{
public:
    int id, uid;
    bool visit;
    std::set<Node *> addjNode;
    std::set<int> addjID, repID, childID;
};

class Graph{
public:
    std::vector<Node*> nodes;
    std::vector<Graph*> connG;
    void printGraph(int mode=0);
};

enum MODE {ACC, ECC};


class clipMa
{
public:
    double **accMa;
    double **eccMa;
    int n,m;

    ~clipMa();
};

typedef std::pair<int,int> pint;

class Database
{
public:
    clip *clips;
    double accThres, eccThres, clipW, clipH;
    double **accMa;
    double **eccMa;
    clipMa **simMa;
    int setN;           // there is setN sampling clip within a marker
    int sN,sM;          // setN is decomposed to sN*sM
    int calCount;
    int samp_nx, samp_ny;

    Graph *graph;
    int N;
    MODE mode;
    std::vector<int> clipSet;
    std::vector<point> clipPos;
    std::vector<std::vector<pint> > clusters;

    Database()
    {
        accThres = 1;
        eccThres = 0;
        clipW = 200;
        clipH = 200;
        mode = ACC;
    }

    void initClips(int Num);

    ~Database();
};

extern Database database;

void _extClip(layout *layoutLib)
{
    printf("extrating clips...\n");
    int clipNum = database.N;
    layoutLib[0].ll.x = (1<<30);
    layoutLib[0].ll.y = (1<<30);
    layoutLib[0].ur.x = -(1<<30);
    layoutLib[0].ur.y = -(1<<30);

    // layoutLib[0].initLayoutllur();
    for (size_t i=0;i<layoutLib[0].polygons.size();i++)
    {
        layoutLib[0].polygons[i]->getLL();
        layoutLib[0].polygons[i]->getUR();
        layoutLib[0].polygons[i]->coordRevise();
        layoutLib[0].polygons[i]->getID(i);
        layoutLib[0].ll.x = std::min(layoutLib[0].polygons[i]->ll.x, layoutLib[0].ll.x);
        layoutLib[0].ll.y = std::min(layoutLib[0].polygons[i]->ll.y, layoutLib[0].ll.y);
        layoutLib[0].ur.x = std::max(layoutLib[0].polygons[i]->ur.x, layoutLib[0].ur.x);
        layoutLib[0].ur.y = std::max(layoutLib[0].polygons[i]->ur.y, layoutLib[0].ur.y);
        // layoutLib[0].polygons[i]->printPolygon();
    }

    layoutLib[0].llYX_sort(); // sort the polygons to enable binary search in chip cutting.
    for (size_t i=0;i<layoutLib[1].polygons.size();i++)
        layoutLib[1].polygons[i]->getCenter();

    layoutLib[1].centerYX_sort(); 

    for (int i = 0; i < clipNum; i++)
    {

        database.clips[i].polygons.clear();
        point ll((1<<30),(1<<30));
        int markerW = 0, markerH = 0;

        for (int j = 0; j < 4; j++)
        {
            ll.x = std::min(ll.x, layoutLib[1].polygons[i]->pt_list[j].x);
            ll.y = std::min(ll.y, layoutLib[1].polygons[i]->pt_list[j].y);
            if(j<4)
            {
                int tw = std::abs(layoutLib[1].polygons[i]->pt_list[j].x - layoutLib[1].polygons[i]->pt_list[j+1].x);
                int th = std::abs(layoutLib[1].polygons[i]->pt_list[j].y - layoutLib[1].polygons[i]->pt_list[j+1].y);
                markerW  = std::max(markerW, tw);
                markerH  = std::max(markerH, th);
            }
        }

        database.clips[i].marker.ll = ll;
        database.clips[i].marker.ur = point(ll.x+markerW,ll.y+markerH);
        database.clips[i].center.x = layoutLib[1].polygons[i]->center.x;
        database.clips[i].center.y = layoutLib[1].polygons[i]->center.y;
        database.clips[i].getBoder();
        extern void clipPolygons_scan(layout *, clip*);
        clipPolygons_scan(layoutLib, &database.clips[i]);
    }
}
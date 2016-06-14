#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_points_boxes(vector <vector <Point*> >hypre_points,int spacing,bool clever,
			  vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
    static int _COUNTER=0;
    int RANK=-1;
    MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
    bool RANKY=RANK==21;
    // if(RANKY)
    cerr << " ENTER BOXES " << RANK << " " << clever << " " << SBoxes.size() << " " << SPoints.size() << endl;
    // SBoxes.clear();
    // SPoints.clear();
    int VOLMIN=4;
    double FILLFACTOR=2.0;
    int MAXY=Misc::pow(2,29);
    int MINY=-Misc::pow(2,29);
    // MPI_Barrier(MPI_COMM_WORLD);
    // int ni=0;
    for(vector <Point*>& hp : hypre_points)
      {
	vector <int> BOX(6);
	vector <int> pos(3);
	for(int B : {0,1,2,3,4,5})
	  BOX[B]=B%2 == 0 ? MAXY:MINY;
	for(Point* &p : hp)
	  {
	    p->get_pos_point(pos);
	    for(int B : {0,2,4})
	      {
		BOX[B]=min(BOX[B],pos[B/2]);
		BOX[B+1]=max(BOX[B+1],pos[B/2]);
	      }
	  }
	Misc::divide(BOX,spacing);
	for(int B : {1,3,5})
	  BOX[B]++;
	KdTree* pHypTree=new KdTree();
	if(clever)
	  pHypTree->LoadKdTree(BOX,hp,spacing,VOLMIN,FILLFACTOR);
	else
	  pHypTree->LoadKdTree(BOX,hp,spacing,-1,2.0);
	int TotalPoints=0;
	int TotalBoxes=0;
	pHypTree->DisplayTree(TotalPoints,TotalBoxes);
	pHypTree->CollectBoxesPoints(SBoxes,SPoints);
	delete pHypTree;
      }
    for(vector<int>& SB : SBoxes)
      {
	for(int B : {1,3,5})
	  SB[B]--;
	Misc::times(SB,spacing);
      }
    _COUNTER++;
    // MPI_Barrier(MPI_COMM_WORLD);
    cerr << " end of boxes " << RANK << " " << _COUNTER << endl;
    cerr << " EXIT BOXES A " << RANK << " " << clever << " " << SBoxes.size() << " " << SPoints.size() << endl;
    if((VOLMIN > 1 || FILLFACTOR < 1.0) && clever)
      any_overlaps(spacing,SBoxes,SPoints);
    // if(RANKY)
    cerr << " EXIT BOXES B " << RANK << " " << clever << " " << SBoxes.size() << " " << SPoints.size() << endl;
  }
}

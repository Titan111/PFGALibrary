#include"PFGA.hpp"

using namespace std;

void init(Genome<int>& g)
{
	int size=g.GetSize();
	for(int i=0;i<size;i++)
		g[i]=i;
}

int main()
{
	Genome<int> g(10);
	init(g);
	cout<<g<<endl;
	return 0;
}


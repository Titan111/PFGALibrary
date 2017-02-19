#include"PFGA.hpp"

using namespace std;

void Cross(Genome<int>* c1,Genome<int>* c2,Genome<int>* p1,Genome<int>* p2)
{
	cerr<<"enter Cross"<<endl;
	int size_=10;
	int cut_n=rand()%size_;

	std::vector<int> x;
	GATool::RandIndex(x,size_);
	cerr<<"pass make x"<<endl;
	std::sort(&x[0],&x[cut_n]);//check
	std::cerr<<"pass sort"<<endl;
	//c1=new Genome<int>(size_);c2=new Genome<int>(size_);
	bool select=true;
	int j=0;

	for(int i=0;i<size_;i++)
	{
		if(i==x[j]&&j<cut_n)
		{
			j++;
			select=!select;
		}

		c1->Gene(i)=select?p1->Gene(i):p2->Gene(i);
		c2->Gene(i)=!select?p1->Gene(i):p2->Gene(i);

	}
	cerr<<c1<<endl;
	cerr<<c2<<endl;
	cerr<<"end cross"<<endl;
}
void init(Genome<int>* g)
{
	int size=g->GetSize();
	for(int i=0;i<size;i++)
		g->Gene(i)=rand()%10;
}

void TenMutation(Genome<int>* g)
{
	int size=g->GetSize();
	int mut_n=rand()%size;
	std::vector<int> x;
	GATool::RandIndex(x,size);
	std::sort(&x[0],&x[mut_n]);//check
	int j=0;
	for(int i=0;i<size;i++)
	{
		if(i==x[j]&&j<mut_n)
		{
			g->Gene(i)=(g->Gene(i)+1)%10;
			j++;
		}
	}
}

double sin_eval(Genome<int>* g)
{
	double ans=0;
	int size=g->GetSize();
	for(int i=0;i<size;i++)
	{
		ans+=sin(i) * g->Gene(i);
	}
	return ans;
}

int main()
{
	/*
	   std::cout<<"function test"<<endl;
	   Genome<int> g1(10);
	   Genome<int> g2(10);
	   Genome<int> *c[2];
	   c[0]=new Genome<int>(10);
	   c[1]=new Genome<int>(10);
	   Cross(c[0],c[1],&g1,&g2);
	   cout<<g1<<endl;
	   cout<<g2<<endl;
	   cout<<*c[0]<<endl;
	   cout<<*c[1]<<endl;
	   cout<<endl;
	 */
	PFGA< Genome<int> > ga;
	ga.size_=10;
	ga.Init=init;
	ga.Mutation=TenMutation;
	ga.Eval=sin_eval;
	int n;
	cin>>n;
	for(int i=0;i<n;i++)
	{
		ga.Step();
		cout<<ga.GetMax()<<endl;
	}
	return 0;
}


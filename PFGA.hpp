#ifndef __PFGA__
#define __PFGA__
#include<algorithm>
#include<iterator>
#include<vector>
#include<iostream>
#include<stdlib.h>
#include<math.h>

using namespace std;

#ifdef DBG
#define ERR_OUT(x) cout<<"\tdbg "<<x<<endl
#define VAR_OUT(x) cout<<"\tdbg "<<#x<<x<<endl
#else
#define ERR_OUT(x)
#define VAR_OUT(x)
#endif

template <class T>
class Genome
{
	private:
		int size_;
		double score_;
		T* data_;
	public:
		Genome<T>()
		{
			score_=0;
		}

		Genome<T>(int s)
		{
			size_=s;
			data_=new T[s];
			score_=0;
		}
		Genome<T>(Genome<T>& g)
		{
			size_=g.GetSize();
			data_=new T[size_];
			score_=0;
			for(int i=0;i<size_;i++)
				data_[i]=g.Gene(i);
		}

		int GetSize(){return size_;}
		void SetScore(double (*eval)(Genome<T>*)){score_=eval(this);}
		double GetScore(){return score_;}

		T& Gene(int n)
		{
			return data_[n];
		}

};

template <typename T>
std::ostream& operator <<(std::ostream& os,Genome<T> g)
{
	int size=g.GetSize();
	os<<"genome:";
	for(int i=0;i<size;i++)
		os<<g.Gene(i)<<" ";
	return os;
}

class GATool
{
	public:
		static void RandIndex(vector<int>& ans,int s)
		{
			for(int i=0;i<s;i++)
				ans.push_back(i);

			for(int i=s-1;i>1;i--)
			{
				int j=rand()%(i);
				int tmp = ans[i];
				ans[i]=ans[j];
				ans[j]=tmp;
			}
		}
};

template <class G>
class PFGA
{
	private:
		std::vector<G*> group_;
		enum MoveCase{A,B,C,D};
	public:
		//setting todo goto private
		int size_;
		int is_fixed_;

		void (*Init)(G*);
		double (*Eval)(G*);
		void (*Mutation)(G*);
		//		void (*Closs)(G*,G*,G*,G*);

		void Cross(G* c1,G* c2,G* p1,G* p2)
		{
			int cut_n=rand()%size_;

			std::vector<int> x;
			GATool::RandIndex(x,size_);
			sort(&x[0],&x[cut_n]);//check
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
		}


		void VecErase(vector<G*>& vec,int n)
		{
			delete vec[n];
			vec.erase(vec.begin()+n);
		}

		int Nnot(int n){return (n+1)%2;}

		G* NewG(int s)
		{
			G* ans;
			ans = new G(s);
			Init(ans);
			return ans;
		}

		int Step()
		{
			ERR_OUT("process prim setting");
			if(group_.size()==0)
			{
				group_.push_back(NewG(size_));
			}

			if(group_.size()==1)
			{
				group_.push_back(NewG(size_));
			}


			ERR_OUT("process cross");
			std::vector<int> p;
			GATool::RandIndex(p,group_.size());
			Genome<int>** child=new Genome<int>*[2];
			child[0]=new G(size_);
			child[1]=new G(size_);
			Cross(child[0],child[1],group_[p[0]],group_[p[1]]);

			ERR_OUT("process mutation");
			int n=rand()%2;
			//cout<<"mutation n:"<<n<<endl;
			//child[n]->SetScore(Eval);
			//cout<<"child:"<<*child[n]<<","<<child[n]->GetScore()<<endl;
			Mutation(child[n]);
			//child[n]->SetScore(Eval);
			//cout<<"child:"<<*child[n]<<","<<child[n]->GetScore()<<endl;

			ERR_OUT("process eval");
			group_[p[0]]->SetScore(Eval);
			group_[p[1]]->SetScore(Eval);
			child[0]->SetScore(Eval);
			child[1]->SetScore(Eval);

			ERR_OUT("process case decision");
			int max_p = (group_[p[0]]->GetScore() > group_[p[1]]->GetScore())?0:1;
			int not_max_p = (max_p+1)%2;
			//TODO rename max and min
			double max=group_[p[max_p]]->GetScore();
			double min=group_[p[not_max_p]]->GetScore();
			int max_c = (child[0]->GetScore() > child[1]->GetScore())?0:1;


			bool out_case=false;
			//TODO make simple
			if(max<=child[0]->GetScore()&&max<=child[1]->GetScore())
			{
				if(out_case)cout<<"case A"<<endl;
				VecErase(group_,p[not_max_p]);
				group_.push_back(child[0]);
				group_.push_back(child[1]);
			}
			else if(min>=child[0]->GetScore()&&min>=child[1]->GetScore())
			{
				if(out_case)cout<<"case B"<<endl;
				VecErase(group_,p[not_max_p]);
				delete child[0];
				delete child[1];
			}
			else if(max>=child[max_c]->GetScore())
			{
				if(out_case)cout<<"case C"<<endl;
				VecErase(group_,p[not_max_p]);
				group_.push_back(child[max_c]);
				delete child[Nnot(max_c)];
			}
			else if(max<=child[max_c]->GetScore())
			{
				if(out_case)cout<<"case D"<<endl;
				VecErase(group_,p[0]);

				int d=0;
				if(p[0]<p[1])d=1;

				VecErase(group_,p[1]-d);
				group_.push_back(child[max_c]);
				delete child[Nnot(max_c)];
				group_.push_back(NewG(size_));
			}else
			{
				cerr<<"error another case"<<endl;
				cerr<<"\tmax"<<max<<endl;
				cerr<<"\tmin"<<min<<endl;
				cerr<<"\tchild_max"<<child[max_c]->GetScore()<<endl;
			}
		}

		int Step(int n)
		{
			for(int i=0;i<n;i++)
				Step();
		}

		double GetMax(){
			int size=group_.size();
			double max_v=0;
			for(int i=0;i<size;i++)
			{
				double value=group_[i]->GetScore();
				if(max_v<=value)
				{
					max_v=value;
				}
			}
			return max_v;
		}

		int GetGroupSize()
		{
			return group_.size();
		}
};
#endif

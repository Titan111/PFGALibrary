#ifndef __PFGA__
#define __PFGA__
#include<vector>
#include<iostream>

template <class T>
class Genome
{
	private:
		int size_;
		T* data_;
	public:
		Genome<T>(int s)
		{
			size_=s;
			data_=new T[s];
		}
		Genome<T>(Genome<T>& g)
		{
			size_=g.GetSize();
			data_=new T[size_];
			for(int i=0;i<size_;i++)
				data_[i]=g[i];
		}

		int GetSize()
		{
			return size_;
		}

		T& operator [](int n)
		{
			return data_[n];
		}
};

template <typename T>
std::ostream& operator <<(std::ostream& os,Genome<T> g)
{
	int size=g.GetSize();
	for(int i=0;i<size;i++)
		os<<g[i]<<" ";
	return os;
}

template <class G>
class PFGA
{
	private:
		std::vector<G> group_;
		std::vector<G> family_;
	public:
		//setting todo goto private
		int size_;
		int is_fixed_;

		int Step();
		int Step(int);
		double (*Eval)(G);
		G (*Closs)(G,G);
		void (*Init)(G&);
		void (*Mutation)(G&);
};
#endif

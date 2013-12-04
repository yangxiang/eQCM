#include <algorithm>
#include <assert.h>
#include <errno.h>
#include <ext/algorithm>
#include <ext/functional>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <math.h>
#include <list>
#include <queue>
#include <string>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <cstdlib>
#include "BitDb.hh"
#include "WG.cpp"
#include "CartesianProductDb.hh"
#include "WGCStatic.hh"


typedef float WeightMeasure;

struct weightcomp {
  bool operator() (const WeightMeasure& lhs, const WeightMeasure& rhs) const
  {return lhs>rhs;}
};

struct sizecomp {
  bool operator() (const Transactionset::size_type& lhs, const Transactionset::size_type& rhs) const
  {return lhs>rhs;}
};

struct edgeinfo{
	Transactionset::transaction_t v_a;
	Itemset::item_t v_b;
};

static void usage() {
	std::cout << "\nUsage:\n"
		"	eQCM [-h] [-g gamma] [-b beta]  [-l lamnda] [-t t_para] input_data_file output_pattern_file\n"
		"Description:\n"
		"	-h	Print the help message.\n"
		"	-g	gamma; default value 0.7; range 0<gamma && gamma<=1.\n"
		"	-b	beta: default value 1; range 0<beta && beta<=1.\n"
		"	-l	lamnda: default value 1; range 1<=lamnda && lamnda<=10.\n"
		"	-t	t_para: default value 1; range 1<=t_para && t_para<=10.\n"
		"	-c	converting_to_absolute: default value 0 (do not convert); range converting_to_absolute={0,1}.\n"
		<< std::endl;
}

int main(int argc, char **argv)
{

  if (argc == 1)
  {
	usage();
	return 1;
  }
  
 
  struct timeval tv;
  gettimeofday(&tv, 0);
  srand48(tv.tv_sec + tv.tv_usec);

  WGCStatic::init_timer();

//  std::ostringstream fi_fname_str;
//  fi_fname_str << WGCStatic::ds_fpath;  
  
  std::string output_file;
  
  //default values begin
  double gamma=0.7; //range 0<gamma && gamma<=1
  double beta=1; //range 0<beta && beta<=1
  int lamnda=1; //range 1<=lamnda && lamnda<=10
  int t_para=1; //range 1<=t_para && t_para<=10
  int converting_to_absolute=0; // range converting_to_absolute={0,1}
  //default values end
  
  int input_para_counter=1;
  while(input_para_counter<argc){
	if(strcmp("-h",argv[input_para_counter])==0){
		usage();
		return 1;
	}
	
	if(strcmp("-g",argv[input_para_counter])==0){
		input_para_counter++;
		gamma=atof(argv[input_para_counter++]);
	}
	else if(strcmp("-b",argv[input_para_counter])==0){
		input_para_counter++;
		beta=atof(argv[input_para_counter++]);
	}
	else if(strcmp("-l",argv[input_para_counter])==0){
		input_para_counter++;
		lamnda=atoi(argv[input_para_counter++]);
	}
	else if(strcmp("-t",argv[input_para_counter])==0){
		input_para_counter++;
		t_para=atoi(argv[input_para_counter++]);
	}
	else if(strcmp("-c",argv[input_para_counter])==0){
		input_para_counter++;
		converting_to_absolute=atoi(argv[input_para_counter++]);
	}	
	else{
		WGCStatic::ds_fpath= argv[input_para_counter++];
		output_file=argv[input_para_counter++];
		break;
	}
  }
  
  //range check begin
  if (!(0<gamma && gamma<=1 && 0<beta && beta<=1 && 1<=lamnda && lamnda<10 && 1<=t_para && t_para<10 && (converting_to_absolute==0 || converting_to_absolute==1) ))
  {
	std::cout<<"At least one input parameter is out of range"<<std::endl;
	usage();
	exit(-1);
  }
  //range check done
  

  
  //managing input and output files
  const std::string fi_fname(WGCStatic::ds_fpath);

  FILE *fi_file;
  if ((fi_file = fopen(fi_fname.c_str(), "rb")) == NULL)
  {
    std::cout << "Dataset " << fi_fname << " does not exist. Exit.\n";
    exit(1);
  }
  
   std::ostringstream out_wcdb;
   out_wcdb<<output_file;
   std::ofstream wcdb_outp(out_wcdb.str().c_str()); 
	if(!wcdb_outp)
	{
		std::cout<<"Fail to create output file. Exit."<<std::endl;
		exit(-1);
	}
	//managing input and output files done

  std::cout << "Skim dataset file to determine its size: " << fi_fname << "\n";
  
  uint64_t Number_of_Rows=0;
  uint64_t Number_of_Cols=0;
  
  size_t line_size = 0;
  char *fi_line = NULL;

  if (getline(&fi_line, &line_size, fi_file)!=-1)
  {
	Number_of_Rows++;
    int item_digit_len;
	WeightMeasure item_weight;
    char *fi_line_scan = fi_line;
	//std::cout<<"reach 0"<<std::endl;
    while (sscanf(fi_line_scan, "%f%n", &item_weight, &item_digit_len) == 1)
	{
		//std::cout<<"item_weight:"<<item_weight<<",";
		//std::cout<<"item_digit_len:"<<item_digit_len<<",";
		fi_line_scan += item_digit_len;
		Number_of_Cols++;
	}
	//std::cout<<"reach 1"<<std::endl;
	while((getline(&fi_line, &line_size, fi_file)!=-1))
		Number_of_Rows++;
	//std::cout<<"reach 2"<<std::endl;
  }

  if (NULL != fi_line)
  {
    free(fi_line);
  }

  if (fclose(fi_file) == -1)
  {
    perror("fclose fi_file");
    exit(1);
  }  
  
  if ((Number_of_Rows<=0)||(Number_of_Cols<=0))
  {
	std::cout<<"The dataset format is incorrect. First Col or Row is not detected. Exit"<<std::endl;
	exit(1);
  }
  WGCStatic::min_item=static_cast<Itemset::item_t>(1);
  WGCStatic::max_item=static_cast<Itemset::item_t>(Number_of_Cols);
  WGCStatic::num_transactions = static_cast<Transactionset::transaction_t>(Number_of_Rows);
  
  //Create graph object
  WG<WeightMeasure> graph(WGCStatic::num_transactions, WGCStatic::min_item, WGCStatic::max_item);
  //Create graph object done
  
  if ((fi_file = fopen(fi_fname.c_str(), "rb")) == NULL)
  {

    std::cout << "Dataset " << fi_fname << " does not exist. Exit.\n";

    exit(1);
  }

  std::cout << "Reading dataset file " << fi_fname << "\n";  

  line_size = 0;
  fi_line = NULL;

  uint64_t row_index=1;
  while (getline(&fi_line, &line_size, fi_file) != -1)
  {

    int item_digit_len;
	WeightMeasure weight;

    char *fi_line_scan = fi_line;
	uint64_t col_index=1;
    while (sscanf(fi_line_scan, "%f%n", &weight, &item_digit_len) == 1)
    {
      if (weight < -std::numeric_limits<WeightMeasure>::max() || weight > std::numeric_limits<WeightMeasure>::max())
      {

		
		std::cout << "weight out of range: " << weight << "\n";
		exit(1);
      }
	  else
	  {
		//std::cout<<weight<<",";
		if (converting_to_absolute==1)
			weight=fabs(weight);//option for abs
		if ((row_index<=WGCStatic::num_transactions)&&(int(col_index)<=int(WGCStatic::max_item)-int(WGCStatic::min_item)+1))
			graph.set_weight(row_index,col_index,weight);
		else
		{
			std::cout << "Dataset format is wrong, possibly not a rectangle. Exit \n";
			exit(1);			
		}
			
	  }
      fi_line_scan += item_digit_len;
	  col_index++;
    }
	//std::cout<<std::endl;
	row_index++;
  }

  if (NULL != fi_line)
  {
    free(fi_line);
  }

  if (fclose(fi_file) == -1)
  {
    perror("fclose fi_file");
    exit(1);
  }
  
  	std::cout<<std::endl;
  
  /*
  for (Transactionset::transaction_t i=1; i<=WGCStatic::num_transactions; i++)
  {
	for (Itemset::item_t j=WGCStatic::min_item; j<=WGCStatic::max_item; j++)
	{
		std::cout<<graph.get_weight(i,j)<<",";
	}
	std::cout<<std::endl;
  }
  */
  
 
  //starting symmetric testing
  if (int(WGCStatic::num_transactions)!=int(WGCStatic::max_item)-int(WGCStatic::min_item)+1)
  {
	std::cout << "This is a bipartite graph. Exit \n";
	exit(1);
  }
  
  for (Transactionset::transaction_t i=1; i<=WGCStatic::num_transactions; i++)
  {
	for (Itemset::item_t j=1; j<=i; j++)
	{
		if (fabs(graph.get_weight(i,j)-graph.get_weight(j,i))>0.001*(fabs(graph.get_weight(i,j))+fabs(graph.get_weight(j,i)))+0.001)
		{
				std::cout << "This graph is not an undirected graph. Exit \n";
				std::cout<<"graph.get_weight(i,j)"<<graph.get_weight(i,j)<<std::endl;
				std::cout<<"graph.get_weight(j,i)"<<graph.get_weight(j,i)<<std::endl;
				exit(1);
		}
		else
		{
			WeightMeasure tmp=graph.get_weight(j,i);
			graph.set_weight(i,j,tmp);
		}
	}
  }
  
  std::cout<<std::endl;
  //std::cout<<graph;  
  
  //ending symmetric testing
  
  
  
  //starting approximation algorithm
  BitDb edgeCover(WGCStatic::num_transactions, WGCStatic::min_item, WGCStatic::max_item);//To set and query whether an edge is covered or not.
  
  std::multimap<WeightMeasure, edgeinfo, weightcomp> edgeRank;//multimap structure ranking edges according to weight from large to small.
  

  
  CartesianProductDb wcdb; //wcdb saves all the dense components
  
  for (Transactionset::transaction_t i=1; i<=WGCStatic::num_transactions; i++)
  {
	for (Itemset::item_t j=WGCStatic::min_item; j<=WGCStatic::max_item; j++)
	{
		if (i==j)
		{
			continue;
		}
		
		edgeinfo tmpinfo;
		tmpinfo.v_a=i;
		tmpinfo.v_b=j;
		edgeRank.insert(std::pair<WeightMeasure, edgeinfo>(graph.get_weight(i,j),tmpinfo));
		edgeCover.set_zero(i,j);//To ensure all zeroes at the beginning
	}
  }
  
  WeightMeasure MAXWEIGHT=edgeRank.begin()->first;
  WeightMeasure Cutoff=gamma*MAXWEIGHT;//The cutoff value of low weight edges that is not worth intializing a new cluster
  
 
  
  
  std::set<Transactionset::transaction_t> TotalSelectedVertices;
  for (std::multimap<WeightMeasure, edgeinfo, weightcomp>::iterator edgeit=edgeRank.begin(); edgeit!=edgeRank.end(); edgeit++)
  {
	//std::cout<<"first edge weight:"<<edgeit->first<<std::endl;
	
	if (edgeit->first<Cutoff)
		break;
	
	if (edgeCover.exists(edgeit->second.v_a, edgeit->second.v_b))
		continue;


	
	//std::cout<<"cutoff value:"<<Cutoff<<std::endl;
	
	//mark the edge as selected (symmetric)
	edgeCover.insert(edgeit->second.v_a, edgeit->second.v_b);
	edgeCover.insert(edgeit->second.v_b, edgeit->second.v_a);

	
	//To be compatible with bipartite case, we still use two sets of vertices but they are equal
	Transactionset V_Trans;
	Itemset V_Item;
	
	V_Trans.push_back(edgeit->second.v_a);
	V_Item.push_back(edgeit->second.v_a);
	V_Trans.push_back(edgeit->second.v_b);
	V_Item.push_back(edgeit->second.v_b);

	std::set<Transactionset::transaction_t> SelectedVertices;	
	SelectedVertices.insert(edgeit->second.v_a);
	SelectedVertices.insert(edgeit->second.v_b);
	
	TotalSelectedVertices.insert(edgeit->second.v_a);
	TotalSelectedVertices.insert(edgeit->second.v_b);
	//std::cout<<"First vertex: "<<edgeit->second.v_a<<"; Second vertex: "<<edgeit->second.v_b<<";";
	
	
	WeightMeasure cluster_weight=edgeit->first;
	
	double density=((double)(2*cluster_weight))/((SelectedVertices.size()-1)*SelectedVertices.size());
	
	while(true)
	{
		Itemset::item_t bestnode=-1;
		WeightMeasure Contribution=-1;
		
		for(Itemset::item_t i=WGCStatic::min_item; i<=WGCStatic::max_item; i++)
		{
			if (SelectedVertices.find(i)!=SelectedVertices.end())
				continue;
				
			WeightMeasure tmp_Weight=0;
			


			for(std::set<Transactionset::transaction_t>::iterator Selected_it=SelectedVertices.begin(); Selected_it!=SelectedVertices.end(); Selected_it++)
			{
				tmp_Weight+=graph.get_weight(*Selected_it, i);
		
			}
			if (tmp_Weight>Contribution)
			{
				Contribution=tmp_Weight;
				bestnode=i;
			}
			
		}

		double c_v_C = ((double)(Contribution))/(SelectedVertices.size());
		double alpha_n=1-1/((double)(2*lamnda*(SelectedVertices.size()+t_para)));
		//std::cout<<"c_v_C: "<<c_v_C<<"; alpha_n: "<<alpha_n<<"; density: "<<density<<"; bestnode: "<<bestnode<<"; Contribution: "<<Contribution<<std::endl;
		
		if(c_v_C >= alpha_n * density)
		{
			V_Trans.push_back(bestnode);
			V_Item.push_back(bestnode);
			for(std::set<Transactionset::transaction_t>::iterator Selected_it=SelectedVertices.begin(); Selected_it!=SelectedVertices.end(); Selected_it++)
			{
				edgeCover.insert(bestnode, *Selected_it);
				edgeCover.insert(*Selected_it, bestnode);
			}
			
			SelectedVertices.insert(bestnode);
			TotalSelectedVertices.insert(bestnode);
			cluster_weight+=Contribution;
		}
		else
		{
			break;
		}
		
		if (SelectedVertices.size()==WGCStatic::num_transactions)
			break;
		
		//Update density!!!
		density=((double)(2*cluster_weight))/((SelectedVertices.size()-1)*SelectedVertices.size());
	}
	
	//SelectedVertices size must be no less than 2.
	CartesianProduct cp(V_Item, V_Trans);
    wcdb.push_back(cp);	
	
  }
  
	if(beta>=1)
		std::cout<<"beta>=1, no merge."<<std::endl;
	else
	{
	  //start merging
		//sort cliques from large to small by a multimap structure
	   std::multimap<Transactionset::size_type, CartesianProduct, sizecomp> CliqueRank;//multimap structure ranking cliques according to their size from large to small.
	   for (CartesianProductDb::iterator w_it=wcdb.begin(); w_it!=wcdb.end(); w_it++)
			CliqueRank.insert(std::pair<Transactionset::size_type, CartesianProduct>(w_it->transactionset.size(), *w_it));
	   
	   //move cliques to a list struture for further merging
	   std::list<CartesianProduct> CPlist;
	   for (std::multimap<Transactionset::size_type, CartesianProduct, sizecomp>::iterator wm_it=CliqueRank.begin(); wm_it!=CliqueRank.end(); wm_it++)
	   {
			CPlist.push_back(wm_it->second);
	   }
	   
	   //intializing s and clear out the multimap
	   CliqueRank.clear();
	   
	   //intializing h
	   std::list<CartesianProduct>::iterator h_iter=CPlist.begin();
	   h_iter++;
		  
	   //intializing j
	   std::list<CartesianProduct>::iterator j_iter=CPlist.begin();
	   
	   std::cout<<"start merging\n";
	   printf("beta: %f\n", beta);
	   
	   while(h_iter!=CPlist.end())
	   {
			//std::cout<<"one mergine cycle\n";
			std::set<uint32_t> mergeclique;
			for (Transactionset::iterator Trans_it=h_iter->transactionset.begin(); Trans_it!=h_iter->transactionset.end(); Trans_it++)
			{
				mergeclique.insert(*Trans_it);
			}
			
			for (Transactionset::iterator Trans_it=j_iter->transactionset.begin(); Trans_it!=j_iter->transactionset.end(); Trans_it++)
			{
				mergeclique.insert(*Trans_it);
			}
			
			int overlap=h_iter->transactionset.size()+j_iter->transactionset.size()-mergeclique.size();
			//std::cout<<"overlap: "<<overlap<<std::endl;
			//std::cout<<"beta*std::min(h_iter->transactionset.size(), j_iter->transactionset.size()): "<<beta*std::min(h_iter->transactionset.size(), j_iter->transactionset.size())<<std::endl;
			
			if (overlap>beta*std::min(h_iter->transactionset.size(), j_iter->transactionset.size()))
			{
				Transactionset Trans_tmp;
				Itemset Item_tmp;
				for(std::set<uint32_t>::iterator mc_it=mergeclique.begin(); mc_it!=mergeclique.end(); mc_it++)
				{
						Trans_tmp.push_back(*mc_it);
						Item_tmp.push_back(*mc_it);
				}
				CartesianProduct cp_tmp(Item_tmp, Trans_tmp);
				CPlist.push_back(cp_tmp);
				mergeclique.clear();
				
				j_iter=CPlist.erase(j_iter);
				
				h_iter=CPlist.erase(h_iter);
				
				j_iter=CPlist.begin();
				
				if(j_iter==h_iter)
					h_iter++;
				
			}
			else
			{
				assert(j_iter!=h_iter);			
				j_iter++;
				
				if(j_iter!=h_iter)
					continue;
				else
				{
					h_iter++;
					j_iter=CPlist.begin();
				}

			}
	   }
	   
	   std::cout<<"end merging\n";
	 
	 
	   wcdb.clear();
	   
	   for(std::list<CartesianProduct>::iterator list_it=CPlist.begin(); list_it!=CPlist.end(); list_it++)
	   {
			wcdb.push_back(*list_it);
	   }
	   
	   CPlist.clear();
	  //end merging
	}
  
  
  //std::cout<<std::endl;
  
  for (CartesianProductDb::iterator w_it=wcdb.begin(); w_it!=wcdb.end(); w_it++)
  {
	for (Transactionset::iterator Trans_it=w_it->transactionset.begin(); Trans_it!=w_it->transactionset.end(); Trans_it++)
	{	
		wcdb_outp<<*Trans_it<<" ";
	}
	//std::cout<<std::endl;	
	
	WeightMeasure cluster_weight=0;
	
	for (Transactionset::iterator Trans_it=w_it->transactionset.begin(); Trans_it!=w_it->transactionset.end(); Trans_it++)
	{	
		for(Itemset::iterator Item_it=w_it->itemset.begin(); Item_it!=w_it->itemset.end(); Item_it++)
		{
			//std::cout<<graph.get_weight(*Trans_it, *Item_it)<<" ";
			if (*Trans_it!=*Item_it)
				cluster_weight+=graph.get_weight(*Trans_it, *Item_it);//each edge weight has been calculated twice
		}
		//std::cout<<std::endl;
	}
	wcdb_outp<<"; Density: "<<(1.0*cluster_weight)/(w_it->transactionset.size()*(w_it->transactionset.size()-1))<<std::endl;
	//std::cout<<std::endl;
  }
   
  wcdb_outp.close();
  //ending approximation algorithm


  return 0;

}

eQCM Program

This program is used for paper "Predicting Glioblastoma Prognosis Networks Using Weighted Gene Co-expression Network Analysis on TCGA Data", Y. Xiang, C.Q. Zhang, K. Huang, BMC Bioinformatics, 2012, 13(Suppl 2), S12. URL: http://www.biomedcentral.com/1471-2105/13/S2/S12/

Software License Agreement: You may use or modify this computer program for research purposes, provided that you properly cite our paper in publication. This computer program is provided on an as is basis and there is no guarantee on the program nor additional support offered. Neither the author(s) nor their institute(s) is liable under any circumstances. This program archive (including this license agreement) may be updated without further notice.

The following of this README file describes how to run the program.

Please update the program to a Linux machine with PBS installed. Type make and obtain an executable file: eQCM.

Sample datasets: 
WeightedUndirectedGraph0.dat
WeightedUndirectedGraph1.dat


Usage: ./bin/eQCM [optional parameters] input_dataset output_file choice

To see parameter hint, simply type ./bin/eQCM -h

Parameters:

input_dataset: The dataset name you want to mine. It should follow the format of the sample dataset. It must be a symmetric matrix.

output_file: The name of the output file. Program will write results into this file. The results are in the form of clusters. One line is one cluster. 
			The number in the cluster is the vertex number. The smallest number is 1. The vertex number corresponds to the line and column of the input matrix,
			i.e., first line and first column of the input matrix correspond to vertex 1, and so on.

		
-g (gamma): This is a parameter of the program (default value 0.7). It should be a float number greater than 0 but no larger than 1. To get some dense components, it is suggested to set between 0.5-1.
		The smaller the value is, the more clusters you get, but the more noise you get too. Suggested value range: 0.7-0.95.
		
-b (beta): This is a parameter of the program (default value 1.0). It should be a float number greater than 0 but no larger than 1. The default value is 1 so that you do not need to merge anything.
		If you would only like to get rid of redundant information, try 0.999999. If you would like to merge more clusters, 0.9 is a suggested value. We would suggest never try a value less than 0.7. Otherwise the results contain too much noise.
		
-l (lamnda): Default value is 1. We would suggest not changing this value. 

-t : Default value 1. We would suggest not changing this value.
		
-c : converting_to_absolute: default value 0 (do not convert). You can specify either 0 or 1, but not other numbers. 0 means original, i.e., each edge weight (each matrix entry) will be used as is in the program. 
		1 means absolution value, i.e., each edge weight (each matrix entry) will be converted to its absolution weight for the program. 
		
		
		
Therefore,a typical command line would look like:
./bin/eQCM ./datasets/WeightedUndirectedGraph0.dat ./results/WeightedUndirectedGraph0_patterns.txt

In order to use the program on clusters, you may also need to write a script file. We have attached two sample script files for you to use.



/*This program uses Rcpp attributes to import pure C++ functions into R. The main
 function takes in an anrmadillo matrix and parses through its rows. NA values are 
 skipped and non-NA values are stored in their respective "group" vectors. Then the 
 coefficient of variation (cv) is taken of each group, using the function calculate_cv.
 Next the function calculate_pool_cv is applied to all the group cv's and produces
 a pooled coefficient of variation for a single row of the input matrix. This 
 process is repeated on all rows of the matrix, the results are stored in a list, 
 and this list is the output of this function.*/

/*The fuction group_size, takes one input, a vector of strings called group. This
 vector is an ordered character vector indicating what "group" (factor level) a specific
 data element belongs to. This function makes a copy of group, called temp. Next a 
 unique function is applied to temp, exposing the levels of the factor. Next we 
 iterate through temp and count how many elements belong to each factor level, these
 counts are stored in gsize.*/

/* The function calculate_cv takes a vector of doubles(group) and calculates the 
 coefficient of variation for that group. This is done by calculating the standard
 deviation of the group and also the mean of the group. Then a ratio of standard 
 deviation to mean is returned.*/

/* The function calculate_pool_cv takes two vectors, the first (cv) stores the 
 coefficient of variation of every group in a row of data. The second vector
 (non_na_values) stores the number of non-missing values per group for one 
 row of data. This function returns the pooled coefficient of variation for a 
 single row of data. The first for loop sums over all groups, the product of 
 coefficient of variation and number of non-missing. The second for loop just 
 adds up all the non-missing values in a row. Then the quotient of these two 
 sums is returned.*/

/* The "main" function is called poolcv, it takes two inputs, a data matrix 
 (armadillo matrix) and a vector of strings which is an ordered character vector 
 indicating what "group" (factor level) a specific data element belongs to.The three 
 nested for loops parse through the matrix of data, row by row. We start at the first 
 row, using the group_size function output, the inner for loops iterate through a 
 row depending on number of groups and size of each group. NA values are skipped and 
 non-NA values are stored in a vector of vectors named "groups." Groups contains the 
 same number of vectors as there is groups per row. After the parsing part of the code, 
 there is another for loop that iterates over "Groups" and calculates the coefficient 
 of variation for each group in a row. If a group contains less than two non-NA values, 
 it's assigned a cv value of zero. Then all the cv's are stored in a vector called groupcv. 
 Next another for loop iterates over "Groups" and counts the number of non-NA values
 per group. If a group has less than two non-NA values, it's assigned a size value of zero. 
 Then all the counts are stored in a vector named nonmiss. Lastly the vectors groupcv 
 and nonmiss are fed into the function calculate_pool_cv which returns the pooled cv for 
 a single row, then pool_cv is stored in a list called "final." This process is applied
 to every row of the matrix. At the end the list "final" of the pooled cv for each row 
 is returned. */


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

std::vector<int> group_size(std::vector<std::string> group)
{
  std::vector<std::string> temp;
  temp = group;
  
  //std::sort(temp.begin(),temp.end());
  
  temp.erase( std::unique( temp.begin(), temp.end() ), temp.end() );
  int tempsize = 0;
  
  tempsize = temp.size();
  
  std::vector<int> gsize(tempsize);
  
  for(unsigned int i = 0;i<group.size();i++)
  {
    for(unsigned int k=0;k<temp.size();k++)
    {
      if(group[i]==temp[k])
        gsize[k]++;
    }
  }
  
  return gsize;
}

double calculate_cv(std::vector<double> numbers)
{
  double mean = 0, sum = 0;
  for (unsigned int i = 0; i < numbers.size(); i++)
  {
    sum = sum + numbers[i];
    mean = sum / numbers.size();
  }
  double var = 0;
  for (unsigned int i = 0; i < numbers.size(); i++)
  {
    var = var + pow((numbers[i] - mean), 2);
  }
  
  return (sqrt(var / (numbers.size() - 1))) / mean;
}

double calculate_pool_cv(std::vector<double> cv, std::vector<double> non_na_values)
{
  double total = 0, sum = 0;
  
  for (unsigned int i = 0; i < cv.size(); i++)
  {
    total = total + cv[i] * non_na_values[i];
  }
  
  for (unsigned int i = 0; i < cv.size(); i++)
  {
    sum = sum + non_na_values[i];
  }
  
  return total / sum;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
std::list<double> pooled_cv_rcpp(arma::mat mtr,std::vector<std::string> group)
{
  std::vector<int> gsize;
  
  gsize = group_size(group); 
  
  std::list<double> final;
  std:: vector<double> cpy(mtr.n_cols);
  std::vector<std::vector<double> >groups;
  
  double tempcv = 0;
  int tempnonmis = 0;
  std::vector<double> groupcv;
  std::vector<double> nonmiss;
  std::vector<double> temp;
  double pool_cv = 0;
  
  for(unsigned int i = 0; i <mtr.n_rows;i++)
  {
    cpy = arma::conv_to< std::vector<double> >::from(mtr.row(i));
    
    for (unsigned int i = 0, j = 0; i < gsize.size(); i++)
    {
      for (unsigned int k = 0; k < gsize[j]; k++)
      {
        if (ISNAN(cpy[k]))
          continue;
        else 
        {
          temp.push_back(cpy[k]);
        }
      }
      groups.push_back(temp);
      cpy.erase(cpy.begin(),cpy.begin()+gsize[j]);
      temp.clear();
      j++;
    }
 
  //calculating cv for each group and pushing_back into gcv vector
  for (unsigned int i = 0; i < groups.size(); i++)
  {
    if(groups[i].size()<= 1)
    {
      groupcv.push_back(NA_INTEGER);
    }
    
    else
    {
      tempcv = calculate_cv(groups[i]);
      groupcv.push_back(tempcv);
    }
  }
  
  //counting the number of nonmissing values per group and pushing_back into nonmiss vector
  
  for (unsigned int i = 0; i < groups.size(); i++)
  { 
       if(groups[i].size()<=1)
        {
           nonmiss.push_back(0);
        }
       else
        {
           tempnonmis = groups[i].size();
           nonmiss.push_back(tempnonmis);
        }
  }
    pool_cv = calculate_pool_cv(groupcv, nonmiss);

    final.push_back(pool_cv);
    
    // clearing these vectors so they will be ready to store the next row's information
    nonmiss.clear();
    groups.clear();
    groupcv.clear();
  }
    
return final;
}
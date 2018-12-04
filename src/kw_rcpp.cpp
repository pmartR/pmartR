/*This program uses Rcpp attributes to import pure C++ functions into R. The main
function takes in an anrmadillo matrix(which has already had the anova
filter applied to it) and parses through its rows. NA values are skipped 
and non-NA values are stored in a vector of vectors called groups, depending
on their group membership. All the non-NA values are concatenated into a vector which
is sorted in ascending order. Each non-NA value is assigned a rank, all the ranks
are stored in another vector called rankvec. The ranks are used in the 
calculation of the Kruskall Wallis h statistic to then calculate the p-value of 
that specific row.*/

/*The fuction group_size takes one input, a vector of strings called group. This
vector is an ordered character vector indicating what "group" factor level a 
specific data element belongs to. This function makes a copy of group, called temp.
Next a unique function is applied to temp, exposing the levels of the factor. 
Next we iterate through temp and count how many elements belong to each level, 
these counts are stored in gsize.*/

/*The function calculate_kwh takes two inputs, a vector of ranks and a vector 
of the number of non-NA values per group.This function implements the 
mathematical formula for Kruskal Wallis test statistic. It returns a numeric 
value(p-value).*/

/*The function compute_pvalue, simply takes in the Kruskall Wallis test statistic
and computes a p-value.*/

/*The main function is the kwh function, it takes two inputs, an armadillo 
matrix and a vector of strings called "group". The vector "group" is fed to
the group_size function which will return a vector of group sizes. One row 
of the armadillo matrix is stored in a vector named "cp". We iterate over
"cp" and collect all non-NA values and store them in a vector of vectors
called "groups". We iterate over "groups" and take the size of each group 
and store it in a vector called "nonmiss".Next we concatenate all the 
vectors from "groups" into one large vector called "cpy." Next we sort 
"cpy" by its indices and store the indices in the vector "yvec." Then
we allocate a vector "rankvec", the same size as "yvec". We assing a rank
to each of the elements of "cpy" and keep them in the same order as "yvec".
We then take "rankvec" and and un-concatenate the vector into three smaller
vectors which are stored in "ranks" (the ranks of each group). Finally "ranks"
is input into calculate_kwh and a p-value is returned for that specific row
of the armadillo matrix. This process is repeated untill the last row of the 
matrix is reached, all the p-values are stored in order in a list named
"final."  */

#include <RcppArmadillo.h>
#include <boost/math/distributions/chi_squared.hpp>
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
using namespace Rcpp;

std::vector<int> gp_size(std::vector<std::string> group)
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

double calculate_kwh(std::vector <std::vector <double> > ranks, std::vector <double> nonmiss_sizes)
{
  std::vector <double> rsums;
  double tempsum = 0,temporary = 0;
  double n = 0,summation = 0, result = 0; 
  
  for(unsigned int i = 0;i<nonmiss_sizes.size();i++)
  {
    n = n + nonmiss_sizes[i];
  }
  
  for(unsigned int i = 0,j=0;i<ranks.size();i++)
  {
    for(int k = 0;k<nonmiss_sizes[j];k++)
    {
      tempsum = tempsum + ranks[i][k];
    }
    
    temporary = pow(tempsum,2)/ nonmiss_sizes[i];
    rsums.push_back(temporary);
    tempsum = 0;
    temporary = 0;
    j++;
  }
  
  for(unsigned int i = 0;i<rsums.size();i++)
  {
    summation = summation + rsums[i];
  }
  
  result = ((12/(n*(n+1)))*summation)- 3*(n+1);
  
  return result;
}

double compute_pvalue(double h,std::vector<double> nonmiss)
{
  double p = 0;
  int df = 0; 
  df = nonmiss.size()-1;
  
  boost::math::chi_squared mydist(df);
  p = boost::math::cdf(mydist, h);
  
  return 1 - p;
}

// [[Rcpp::export]]
std::list<double> kw_rcpp(arma::mat mtr,std::vector<std::string> group)
{
  std::vector<std::string> unique_groups = group; 
  std::sort(unique_groups.begin(),unique_groups.end());
  unique_groups.erase( std::unique( unique_groups.begin(), unique_groups.end() ), unique_groups.end() );
  
  std::list<double> final;
  std:: vector<double> cp(mtr.n_cols);
  std::vector<std::vector<double> >groups;
  
  double tempnonmis = 0;
  
  std::vector<double> nonmiss;
  std::vector<double> tmp;
  
  for(unsigned int i = 0; i <mtr.n_rows;i++)
  {
        cp = arma::conv_to< std::vector<double> >::from(mtr.row(i));
        
        for (unsigned int i = 0; i < unique_groups.size(); i++)
        {
          for (unsigned int k = 0; k < cp.size(); k++)
          {
            if (ISNAN(cp[k])){
              continue;
            }
            else if(unique_groups[i] == group[k]) {
              tmp.push_back(cp[k]);
            }
          }
          
          groups.push_back(tmp);
          tmp.clear();
        
       }
    
            int all_na_count = 0;
    
            for (unsigned int i = 0; i < groups.size(); i++)
            {
              if(groups[i].size()==0)
              {
                all_na_count++;
              }
    
             else
             {
               tempnonmis = groups[i].size();
               nonmiss.push_back(tempnonmis);
             }
           }

            if((groups.size()-all_na_count) < 2)
            {
              final.push_back(NA_REAL);
              nonmiss.clear();
              groups.clear();
            }
    
          else
            {
              std::vector<double> cpy;
      
              for(unsigned int i = 0; i <groups.size();i++)
              {
                if(groups[i].size() == 0)
                {
                  continue;
                }
                
                else
                {
                  cpy.insert(cpy.end(), groups[i].begin(), groups[i].end());
                }
              }
      
              std::vector <double> rankvec (cpy.size());
              std::vector<double> yvec(cpy.size());
              std::size_t n(0);
      
              std::generate(yvec.begin(), yvec.end(), [&]{ return n++; });
              std::sort(yvec.begin(), yvec.end(), [&](int i1, int i2) {return cpy[i1] < cpy[i2]; } );
      
              for(unsigned int i = 0,j = 1;i<yvec.size();i++)
              {
                rankvec[yvec[i]]= j;
                j++;
              }
      
              double tempr = 0;
              std::vector<double>tempvec;
              std::vector<std::vector<double> > ranks;
      
                for (unsigned int i = 0, j = 0; i <nonmiss.size(); i++)
                {
                    for (unsigned int i = 0; i < nonmiss[j]; i++)
                    {
                      tempr = rankvec[i];
                      tempvec.push_back(tempr);
                    }
        
                  ranks.push_back(tempvec);
                  rankvec.erase(rankvec.begin(),rankvec.begin()+nonmiss[j]);
                  tempvec.clear();
        
                  j++;
                }
      
                double h = 0, p = 0;
      
                h = calculate_kwh(ranks,nonmiss);
                p = compute_pvalue(h,nonmiss);
                
                final.push_back(p);
      
                nonmiss.clear();
                groups.clear();
                ranks.clear();
                rankvec.clear();
                yvec.clear();
                cpy.clear();
            }
      }
  
  return final;
}

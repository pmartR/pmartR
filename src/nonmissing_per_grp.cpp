#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


std::vector<int> grp_siz(std::vector<std::string> group)
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

// [[Rcpp::export]]

arma::Mat<int> nonmissing_per_grp(arma::mat mtr,std::vector<std::string> group) 
{
  std::vector<int> nvec;
  nvec = grp_siz(group);
  
  std::vector<double> cpy;
  std::vector<int> non_missing;
  int count = 0;
  
  arma::Mat<int> final(mtr.n_rows,nvec.size());
  arma::Row<int> row;
  
  std::vector<int> temp;
  
  for(unsigned int z = 0; z <mtr.n_rows;z++)
  {
    cpy = arma::conv_to< std::vector<double> >::from(mtr.row(z));
    
    for (unsigned int i = 0, j = 0; i < nvec.size(); i++)
    {
      
      for (unsigned int k = 0; k < nvec[j]; k++)
      {
        
        if(ISNAN(cpy[k]))
          continue;
        
        else count++;
        
      }
      
      non_missing.push_back(count);
      cpy.erase(cpy.begin(),cpy.begin()+nvec[j]);
      count = 0;
      j++;
      
    }
    
    row = arma::conv_to<arma::Row<int> >::from(non_missing);
    
    final.row(z)=row;
  
    
    non_missing.clear();
    row.clear();
    cpy.clear();
    
  } 
  
  return final;
}

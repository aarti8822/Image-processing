#include<iostream>
#include<cmath>
#include"project2.hpp"
using namespace std;
int main(){
int n;
cout<<"Enter the number of rows and columns"<<endl;
cin>>n;
double *a,*b,*c,*v,*s,*d,*u,*f,*g,*h;
a=new double[n*n];
b=new double[n*n];
c=new double[n*n];
v=new double[n*n];
s=new double[n*n];
d=new double[n*n];
u=new double[n*n];
f=new double[n*n];
g=new double[n*n];
h=new double[n*n];
for(int i=0;i<n;i++){
for(int j=0;j<n;j++){
cin>>a[i*n+j];
}
}
cout<<"The given matrix is-"<<endl;
print(a,n,n);
cout<<"\n"<<endl;
transpose(a,b,n,n);
cout<<"The transpose of given matrix is i.e A^T="<<endl;
print(b,n,n);
cout<<"\n"<<endl;
mult(b,a,c,n,n,n);
cout<<"The A^T*A will be ="<<endl;
print(c,n,n);
cout<<"\n"<<endl;
eigenvector(c,v,s,n);
cout<<"The matrix V is-"<<endl;
print(v,n,n);
cout<<"\n"<<endl;
cout<<"The matrix sigma is-"<<endl;
print(s,n,n);
cout<<"\n"<<endl;
mult(a,b,d,n,n,n);
cout<<"The A*A^T will be ="<<endl;
print(d,n,n);
cout<<"\n"<<endl;
eigenvectoru(d,u,n);
cout<<"The matrix U is-"<<endl;
print(u,n,n);
cout<<"\n"<<endl;
transpose(v,f,n,n);
mult(s,f,g,n,n,n);
mult(u,g,h,n,n,n);
cout<<"Multiplying U,Sigma and V^T and we will get given matrix A"<<endl;
print(h,n,n);
cout<<"\n"<<endl;
delete []a;
delete []b;
delete []c;
delete []v;
delete []s;
delete []u;
delete []f;
delete []g;
delete []h;

return 0;
}
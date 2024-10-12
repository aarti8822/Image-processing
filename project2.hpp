// SINGULAR VALUE DECOMPOSITION ETHOD WAS USED IN THIS PROJECT TO REDUCE THE DIMENSION OF IMAGE MATRIX


#include<iostream>
#include<cmath>
using namespace std;


//For printining any a=m*n matrix


void print(double *a,int m,int n){
for(int i=0;i<m;i++){
cout<<"[";
for(int j=0;j<n;j++){
cout<<a[i*n+j]<<"\t";
}
cout<<"]"<<endl;
}
}

//For transpose of any given matrix.If a=m*n matrix,then b is the transpose matrix of and of order n*m.


void transpose(double *a,double *b,int n,int m){
for(int i=0;i<n;i++){
for(int j=0;j<m;j++){
b[i*m+j]=a[j*n+i];
}
}
}

//For matrix multiplication.Let a=m*n,b=n*p matrix then c=m*p matrix


void mult(double *a,double *b,double *c,int m,int n,int p){
for(int i=0;i<m;i++){
for(int k=0;k<p;k++){
c[i*p+k]=0;
for(int j=0;j<n;j++){
c[i*p+k]+=a[i*n+j] * b[j*p+k];
}
}
}
}


// For printing any a=n*1 vector


void vectorp(double *a,int n){
for(int i=0;i<n;i++){
cout<<"[";
cout<<*(a+i);
cout<<"]"<<endl;
}
}

//  For extraction of pth column vector 'b' from given matrix 'a' of m*n order

  
void vector(double *a,double *b,int m,int n,int p){
int s=0;
for(int i=0;i<m;i++){
b[s++]=a[i*n+p];
}
}

// For norm of any b=m*1 vector.


void norm(double *b,int m){
double sum=0;
for(int i=0;i<m;i++){
sum+=b[i] * b[i];
}
sum=sqrt(sum);
for(int i=0;i<m;i++){
b[i]=b[i]/sum;
}
}

//For substraction of any two n*1 vectors


void vectors(double *a,double *b,double *c,int n){
for(int i=0;i<n;i++){
c[i]=a[i] - b[i];
}
}

//For dot product of any two 'n' vectors


double dot(double *a,double *b,int n){
double d=0;
for(int i=0;i<n;i++){
d+=a[i] * b[i];
}
return d;
}

//For multiplying any a=n*1 vector with any scalar 'd'
void scalar(double *a,int n,double d){
for(int i=0;i<n;i++){
a[i]=d * a[i];
}
}

//To find absolute maximum from any n*1 vector


double maxm(double *a,int n){
double max=fabs(a[0]);
for(int i=0;i<n;i++){
if(max<=fabs(a[i])){
max=fabs(a[i]);
}
}
return max;
}


// For dividing any a=n*1 vector with 'z'
void reqv(double *a,int n,double z){
for(int i=0;i<n;i++){
a[i]=a[i]/z;
}
}


// For finding maximum or minimum eigen value and its corresponding eigen vector.Here 'a' is any given matrix of order n*n and 'b' is initial n*1 vector
 

void power(double *a,double *b,int n){
double *c;
c=new double[n*n];
double z=0.0;
double d=0.0;
double e=0.001;
while(fabs(e)>=0.001){
mult(a,b,c,n,n,1);
z=maxm(c,n);
reqv(c,n,z);
for(int i=0;i<n;i++){
b[i]=c[i];
}
mult(a,b,c,n,n,1);
d=maxm(c,n);
reqv(c,n,d);
for(int i=0;i<n;i++){
b[i]=c[i];
}
 e=d-z;
}
cout<<"Eigen value is "<<d<<endl;
cout<<"\n"<<endl;
cout<<"Its corresponding eigen vector"<<endl;
vectorp(b,n);
}

// For finding Matrix 'q' of given matrix 'a' of order m*n.    i.e A=QR;here we will find Q of given matrix A;


void qdecomposition(double *a,double *q,int m,int n){
double *b,*c,*f;
b=new double[m];
c=new double[m];
f=new double[m];
vector(a,b,m,n,0);//  Took out zeroth column from matrix 'a'
norm(b,m);//   Normalized it.
for(int i=0;i<m;i++){// Made the vector b as zeroth column of matrix 'q'.
for(int j=0;j<n;j++){
if(j==0){
q[i*n+j]=b[i];
}
else{
q[i*n+j]=0;
}
}
}

// Here we are trying to find other remaining columns of matrix 'q'.
for(int i=1;i<n;i++){
double *t;
t=new double[m]();
for(int j=0;j<i;j++){
vector(a,b,m,n,i);
vector(q,c,m,n,j);
double z=dot(b,c,m);
scalar(c,m,z);
for(int k=0;k<m;k++){
t[k]+=c[k];
}
}
vector(a,b,m,n,i);
vectors(b,t,f,m);
norm(f,m);
for(int k=0;k<m;k++){
q[k*n+i]+=f[k];
}
}
delete[] b;
delete[] c;
delete[] f;
}

// Here we are finding matrix 'r' of order n*n matrix   i.e  A=QR;here we are Finding R;


void rdecomposition(double *a,double *g, double *h,int m,int n){
double *y;
y=new double[n*m];
 transpose(g,y,n,m);
mult(y,a,h,n,m,n);
delete[] y;
}


// Here we are trying to multply any two matrix and changing of one given matrix into resultant matrix.  i.e  a=a*b;here a and b are matrix;


void nmatrixmul(double *a,double *b,int n){
double *c;
c=new double[n*n];
for(int i=0;i<n;i++){
for(int k=0;k<n;k++){
c[i*n+k]=0;
for(int j=0;j<n;j++){
c[i*n+k]+=a[i*n+j] * b[j*n+k];
}
}
}
for(int i=0;i<n;i++){
for(int j=0;j<n;j++){
a[i*n+j]=c[i*n+j];
}
}
delete[] c;
}

//  Here we are trying to find eigen values and eigen vectors of any given n*n matrix 'a'i.e matrix 'V' and Sigma using qr decomposition.


void eigenvector(double *a,double *x,double *y,int n){
double *u,*b,*q,*r,*s;
u=new double[n*n];
b=new double[n*n];
q=new double[n*n];
r=new double[n*n];
s=new double[n];
for(int i=0;i<n;i++){
for(int j=0;j<n;j++){
if(i==j){
u[i*n+j]=1;
}
else{
u[i*n+j]=0;
}
}
}
for(long int i=0;i<100000;i++){
qdecomposition(a,q,n,n);
nmatrixmul(u,q,n);
rdecomposition(a,q,r,n,n);
mult(r,q,b,n,n,n);
for(int j=0;j<n;j++){
for(int k=0;k<n;k++){
a[j*n+k]=b[j*n+k];
}
}

i++;
}
for(int i=0;i<n;i++){
for(int j=0;j<n;j++){
if(i==j){
*(y+i*n+j)=sqrt(*(a+i*n+j));
}
else{
*(y+i*n+j)=0;
}
}
}
for(int k=0;k<n;k++){
vector(u,s,n,n,k);
norm(s,n);
for(int i=0;i<n;i++){
for(int j=0;j<n;j++){
if(j==k){
*(x+i*n+j)=*(s+i);
}
}
}
}
delete[] u;
delete[] b;
delete[] q;
delete[] r;
delete[] s;
}

// Here we are trying to find matrix 'U' by the same method.
void eigenvectoru(double *a,double*x,int n){
double *u,*b,*q,*r,*s;
u=new double[n*n];
b=new double[n*n];
q=new double[n*n];
r=new double[n*n];
s=new double[n];
for(int i=0;i<n;i++){
for(int j=0;j<n;j++){
if(i==j){
u[i*n+j]=1;
}
else{
u[i*n+j]=0;
}
}
}
for(long int i=0;i<100000;i++){
qdecomposition(a,q,n,n);
nmatrixmul(u,q,n);
rdecomposition(a,q,r,n,n);
mult(r,q,b,n,n,n);
for(int j=0;j<n;j++){
for(int k=0;k<n;k++){
a[j*n+k]=b[j*n+k];
}
}

i++;
}
for(int k=0;k<n;k++){
vector(u,s,n,n,k);
norm(s,n);
for(int i=0;i<n;i++){
for(int j=0;j<n;j++){
if(j==k){
*(x+i*n+j)=*(s+i);
}
}
}
}
delete[] u;
delete[] b;
delete[] q;
delete[] r;
delete[] s;
}








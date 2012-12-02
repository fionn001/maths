#include <iostream>
#include <iomanip>  
#include <stdlib.h> 
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
using namespace std;
void shuffle(vector<int>& shuff, vector<int>& li, int x);
void list(vector<int>& li,int x);
void conf(vector<int>& shuff,vector<int>& li,vector<int>& pos,vector<int>& row,vector<int>& col,int x,int y, int z);
void world(vector<int>& col,vector<int>& row,int x,int y,int z,vector<double>& m,double HOP,double POT);
int pauli_l(vector<int>& V,int rand,int x,int y,int z);
int pauli_r(vector<int>& V,int rand,int x,int y,int z);
int metrop(double init, double fin);
void elements(vector<double>& m,double DT,double HOP,double POT);
double occ_lr(vector<int>& A,int p, vector<double>& m,int x,int y,int z);
double occ_rl(vector<int>& A,int p, vector<double>& m,int x,int y,int z);
double diag(vector<int>& B, int pl, vector<double>&m);
double energy(vector<int>& A,vector<int>& B,int x,int y,int z,vector<double>& m, double HOP,double POT);
int main(){
	int N,T,site,i,j;
	double dt,hop,pot;
	dt=0.1;				//imaginary time
	hop=1;				//hopping parameter
	pot=0;				//V parameter
	N=20;				//number of electrons
	T=21;				//number of time steps
	site=40;				//number of sites
	srand(time(NULL));
	vector<int> row(N*T),col(N*T),pos(N),shuff(site),V(N*T),li;
	vector<double> m(4);
	elements(m,dt,hop,pot);
	for(i=0;i<3;i++){
		cout << m[i] << endl;
	}
	conf(shuff,li,pos,row,col,N,T,site);	
	sort(col.begin(),col.end());	
	world(col,row,site,N,T,m,hop,pot);
}

void world(vector<int>& col,vector<int>& row,int x,int y,int z,vector<double>& m,double HOP,double POT){
	int i,r,j,count1,count2,count3,k;
	double en_tot;
	en_tot=0.0;
	count3=0;
	vector< vector<int> > MAT(z,vector<int>(x,0));
	ofstream file;
	ofstream file2;
	file.open("col1.txt");
	file2.open("row2.txt");
	vector<int> tmp_col(y*z);
	for(j=0;j<100000;j++){
		r=rand()%(z*y);
		if((row[r]==0)||(row[r+1]==(z-1))){
				col[r]=col[r];
				col[r+1]=col[r+1];
		}
		else{
		if(((((col[r]%2)==1)&&((row[r]%2)==0))||(((col[r]%2)==0)&&((row[r]%2)==1)))&&(col[r]==col[r+1])){
			if(pauli_l(col,r,x,y,z)>0){				
				col[r]=col[r];
				col[r+1]=col[r+1];
			}
			else{							//shift world-line left
				tmp_col=col;
				if((tmp_col[r]-1)<0){		
					tmp_col[r]=x-1;
					tmp_col[r+1]=x-1;
					//cout << "left" << col[r] <<"  "<< metrop(occ_rl(col,r,m,x,y,z),occ_rl(tmp_col,r,m,x,y,z)) << endl;
					if(metrop(occ_rl(col,r,m,x,y,z),occ_rl(tmp_col,r,m,x,y,z))>0){
						col=tmp_col;
					}
					else{
						col=col;
					}
				}
				else{
					tmp_col[r]=tmp_col[r]-1;
					tmp_col[r+1]=tmp_col[r+1]-1;
					//cout << "left" << col[r] <<"  "<< metrop(occ_rl(col,r,m,x,y,z),occ_rl(tmp_col,r,m,x,y,z)) << endl;
					if(metrop(occ_rl(col,r,m,x,y,z),occ_rl(tmp_col,r,m,x,y,z))>0){
						col=tmp_col;
					}
					else{
						col=col;
					}
				}
			}
		}
		else if((((col[r]%2)==0)&&((row[r]%2)==0)||((col[r]%2)==1)&&((row[r]%2)==1))&&(col[r]==col[r+1])){
			if(pauli_r(col,r,x,y,z)>0){
				col[r]=col[r];
				col[r+1]=col[r+1];
			}
			else{
				tmp_col=col;				
				tmp_col[r]=(tmp_col[r]+1)%x;
				tmp_col[r+1]=(tmp_col[r+1]+1)%x;
				//cout << "right" << col[r] <<"  "<< metrop(occ_lr(col,r,m,x,y,z),occ_lr(tmp_col,r,m,x,y,z)) << endl;
				if(metrop(occ_lr(col,r,m,x,y,z),occ_lr(tmp_col,r,m,x,y,z))>0){
					col=tmp_col;
				}
				else{
					col=col;
				}
			}
		}
		else{
			col[r]=col[r];
			col[r+1]=col[r+1];
		}
		
	}
	if((j>50000)&&((j%5)==0)){
		en_tot=en_tot+energy(row,col,x,y,z,m,HOP,POT);
	}
	}
	//en_tot=en_tot+energy(row,col,x,y,z,m,HOP,POT);
	cout <<"energy:  "<<en_tot/(x*z*(10000)) << endl;
	/*
	for(i=0;i<(y*z);i++){
			file << col[i] << endl;
			file2 << row[i] << endl;
	}
	for(i=0;i<z;i++){
	for(k=0;k<x;k++){
		MAT[i][k]=0;
	}
	}
	for(i=0;i<z*y;i++){
		MAT[row[i]][col[i]]=1;
	}
	for(i=0;i<z;i++){
	for(k=0;k<x;k++){
		cout<<MAT[i][k]<<" ";
	}
		cout<<endl;
	}
	*/
	cout << endl; 
	cout << endl;
	cout << j << endl;
	file.close();
	file2.close();	
}
int metrop(double init,double fin){
	double prob,r_num;
	prob=fin/init;
	if(prob>10){
	//cout <<" ratio:   " <<prob << endl;
	}
	if(prob>1.0){
		return(1);
	}
	else{
		r_num=((double)rand()/((double)RAND_MAX));
		if(prob>r_num){
			return(1);
		}
		else{
			return(0);
		}
	}
}
double energy(vector<int>& A,vector<int>& B,int x,int y,int z,vector<double>& m, double HOP,double POT){
	int i,j,k;
	double en;
	en=0.0;
	vector< vector<int> > wl(z,vector<int>(x,0));
	for(i=0;i<z;i++){
	for(j=0;j<x;j++){
		wl[i][j]=0;
	}
	}
	for(i=0;i<z*y;i++){
		wl[A[i]][B[i]]=1;
	}
	
	for(i=0;i<(z-1);i=i+2){
	for(j=1;j<(x-1);j=j+2){
		if(wl[i][j]==wl[i][j+1]){
			en=en+POT/4.0;
		}
		else if(wl[i][j]==wl[i+1][j]){
			en=en-(1.0*HOP*(m[1]/m[0])+POT/4.0);
		}
		else{
			en=en-(1.0*HOP*(m[0]/m[1])+POT/4.0);	
		}	
	}
	}
	for(i=1;i<(z-1);i=i+2){
	for(j=0;j<(x-1);j=j+2){
		if(wl[i][j]==wl[i][j+1]){
			en=en+POT/4.0;
		}
		else if(wl[i][j]==wl[i+1][j]){
			en=en-(1.0*HOP*(m[1]/m[0])+POT/4.0);
		}
		else{
			en=en-(1.0*HOP*(m[0]/m[1])+POT/4.0);	
		}	
	}
	}
	return(en);
}
double diag(vector<int>& B, int pl, vector<double>&m){
	int count,count1;
	count=0;
	count1=0;
	if(abs(B[pl-1]-B[pl])==0){				//straight line
		count++;
	}
	else{							//diagonal
		count1++;
	}
	if(abs(B[pl+2]-B[pl])==0){				//straight
		count++;
	}
	else{							//diag
		count1++;
	}
	if(count==2){						//two straigh edges
		return(m[0]*m[0]);
	}
	else if(count1==2){					//two diagonals
		return(m[1]*m[1]);
	}
	else{							//one straight and one diagonal
		return(m[1]*m[0]);
	}
}
double occ_lr(vector<int>& A,int p, vector<double>& m,int x,int y,int z){
	int count1,count2;
	count1=0;
	count2=0;
	//looking at square at i-1
	if((A[(p/z)*z]==A[0])){				//left end point
		if((A[p]-1)<0){
			if((A[p+(y-1)*z]==(x-1))){
				count1++;
			}	
		}
		else if(A[p+(y-1)*z]==(A[p]-1)){
			count1++;
		}
		else{
			count1=count1;
		}
	}
	else{				//right end point
		if((A[p]-1)<0){	
			if((A[p-z])==(x-1)){
				count1++;
			}
		}
		else if(A[p-z]==(A[p]-1)){
			count1++;
		}
		else{
			count1=count1;
		}
	}
	//looking at square at i+2
	if((A[(p/z)*z]==A[(y-1)*z])){				//right end point case
		if(A[p-(y-1)*z]==((A[p]+2)%x)){
			count2++;
		}
	}
	else{
		if(A[p+z]==((A[p]+2)%x)){			//everywhere else
			count2++;
		}
	}
	//matrix elements
	if((count1>0)&&(count2>0)){
		return(m[0]*m[2]*diag(A,p,m));
	}
	else if(((count2==1)&&(count1==0))){
		return(m[0]*m[0]*diag(A,p,m));
	}
	else if(((count2==0)&&(count1==1))){
		return(m[3]*m[2]*diag(A,p,m));
	}
	else{
		return(m[0]*m[3]*diag(A,p,m));
	}
}
double occ_rl(vector<int>& A,int p, vector<double>& m,int x,int y,int z){
	int count1,count2;
	count1=0;
	count2=0;
	//looking at square at i+1
	if((A[(p/z)*z]==A[(y-1)*z])){				//right end point case
		if(A[p-(y-1)*z]==((A[p]+1)%x)){
			count2++;
		}
	}
	else{
		if(A[p+z]==((A[p]+1)%x)){			//everywhere else
			count2++;
		}
	}
	//looking at square at i-2
	if((A[(p/z)*z]==A[0])){				//left end point
		if((A[p]-2)<0){
			if((A[p+(y-1)*z]==(x-2))){
				count1++;
			}	
		}
		else if(A[p+(y-1)*z]==(A[p]-2)){
			count1++;
		}
		else{
			count1=count1;
		}
	}
	else{				//right end point
		if((A[p]-2)<0){	
			if((A[p-z])==(x-2)){
				count1++;
			}
		}
		else if(A[p-z]==(A[p]-2)){
			count1++;
		}
		else{
			count1=count1;
		}
	}
	if((count1>0)&&(count2>0)){
		return(m[0]*m[2]*diag(A,p,m));
	}
	else if(((count2==1)&&(count1==0))){
		return(m[3]*m[2]*diag(A,p,m));
	}
	else if(((count2==0)&&(count1==1))){
		return(m[0]*m[0]*diag(A,p,m));
	}
	else{
		return(m[0]*m[3]*diag(A,p,m));
	}
}
void elements(vector<double>& m,double DT,double HOP,double POT){
	m[0]=cosh(DT*HOP)*exp((DT*POT)/4.0);		//straight
	m[1]=sinh(DT*HOP)*exp((DT*POT)/4.0);		//diag
	m[2]=exp(-1*(DT*POT)/4);			//zero
	m[3]=exp(-1*(DT*POT)/4);			//double
}
int pauli_l(vector<int>& V,int rand,int x,int y,int z){
	int count1;
	count1=0;
	//right endpoint
	if((V[(rand/z)*z]==V[(y-1)*z])){
		if((V[rand]-1)<0){	
			if(((V[rand-z+1])==(x-1))||((V[rand-z])==(x-1))){
				count1++;
			}
		}
		else if(((V[rand+1]-1)==V[rand-z+1])||((V[rand]-1)==V[rand-z])){		
			count1++;
		}	
		else{
			count1=count1;
		}
	}
			//left end point
	else if((V[(rand/z)*z]==V[0])){
		if((V[rand]-1)<0){
			if(((V[rand+(y-1)*z+1])==(x-1))||((V[rand+(y-1)*z])==(x-1))){
				count1++;
			}
		}
		else if(((V[rand+1]-1)==V[rand+(y-1)*z+1])||((V[rand]-1)==V[rand+(y-1)*z])){		
			count1++;
		}	
		else{
			count1=count1;
		}
	}
			//everywhere else
	else{
		if((V[rand]-1)<0){
			if((V[rand-z]==(x-1))||(V[rand-z+1]==(x-1))){
				count1++;
			}
		}
		else if(((V[rand+1]-1)==V[rand-z+1])||((V[rand]-1)==V[rand-z])){		
			count1++;
		}	
		else{
			count1=count1;
		}
	}
	return(count1);
}
int pauli_r(vector<int>& V,int rand,int x,int y, int z){
	int count2;
	count2=0;
	if(V[(rand/z)*z]==V[(y-1)*z]){
		if((((V[rand+1]+1)%x)==V[rand-(y-1)*z+1])||(((V[rand]+1)%x)==V[rand-(y-1)*z])){
			count2++;
		}
		else{
			count2=count2;
		}
	}
	else{
		if((((V[rand+1]+1)%x)==V[rand+z+1])||(((V[rand]+1)%x)==V[rand+z])){
			count2++;
		}
		else{
			count2=count2;
		}
	}
	return(count2);
}
//row and column indices of initial straight line configurations
void conf(vector<int>& shuff,vector<int>& li,vector<int>& pos,vector<int>& row,vector<int>& col,int x,int y, int z){
	int i;
	list(li,z);
	shuffle(shuff,li,z);
	ofstream file;
	ofstream file2;
	//file.open("col.txt",ios::out|ios::app);
	//file2.open("row.txt",ios::out|ios::app);
	for(i=0;i<x;i++){
		pos[i]=shuff[i];
		//cout << pos[i]<<endl;
	}
	for(i=0;i<x*y;i++){
		row[i]=i%y;
		col[i]=pos[i/y];	
		//file<<col[i] << endl;
		//file2<<row[i] << endl;
	}
	for(i=0;i<x*y;i++){
		//cout << col[i] <<"	"<<row[i]<< endl;
	}
	//file.close();
	//file2.close();
}
//randomly assigning initial positions of electrons on lattice
void shuffle(vector<int>& shuff,vector<int> & li, int x){
	int r,i,j;
	vector<int> te;
	
	for(i=0;i<x;i++){
		te.push_back(i);
	}
	j=0;
	for(i=(x-1);i>-1;i--){
		r=(int)(rand()%(i+1));
		shuff[j]=te[r];
		te.erase (te.begin()+r);
		j++;
	}
}
void list(vector<int>& li,int x){
	int i;
	for(i=0;i<x;i++){
		li.push_back(i);
	}
}


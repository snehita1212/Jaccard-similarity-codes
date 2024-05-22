#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/time.h>
#include <omp.h>
#include <pthread.h>
#include <sys/ipc.h>
#include <sys/shm.h>
//#include "shared.h"

// #include <chrono>
using namespace std;
// using namespace std::chrono;

vector<unsigned int> r;
//unordered_map<long long int,vector<long long int>> g_hash;
//unordered_map<long long int, long long int> g_pos;

template<typename T>
vector<T> fileToVec(ifstream& inp){
	int size=0;
	inp>>size;
	vector<T> vec(size);
	unsigned int ele=0;
	for (int i=0; i<size; i++){
		inp>>ele;
		vec[i]=ele;
	}
	//for (auto ele: vec) inp>>ele;
	return vec;
}

vector<unsigned int> pre_hash(string s){

    //kgrams
    vector<string> kgram;
    for (unsigned int i=0; i<s.length(); i++){
        kgram.push_back(s.substr(i,16));
    }
    //kgram to binary
    map<char, string> letter = {{'A', "00"}, {'C', "01"}, {'G', "10"}, {'T', "11"}};
    for (unsigned int i=0; i<kgram.size(); i++) {
        string bin;
        for (char c : kgram[i]) {
            bin += letter[c];
        }
        kgram[i]=bin;
    }
    //binary to int
    vector<unsigned int> hash;
    for (auto& bin: kgram){
        unsigned int num=0;
        for (char c:bin){
            num=num*2+(c-'0');
        }
        hash.push_back(num);
    }
    
    return hash;
}

void win_read(vector<unsigned int> hash, unsigned int k){
    deque<unsigned int> dq;
    // vector<long long int> res;
    for (unsigned int i=0; i<k; i++){
        while(!dq.empty() && hash[i]<=hash[dq.back()]) dq.pop_back();
        dq.push_back(i);
    }
    r.push_back(hash[dq.back()]);
    
    unsigned int ind=dq.front();
    for (unsigned int i=k; i<hash.size(); i++){
        if (dq.front()==i-k) dq.pop_front();
        while(!dq.empty() && hash[i]<=hash[dq.back()]) dq.pop_back();
        dq.push_back(i);
        if (ind!=dq[0]){
            r.push_back(hash[dq.front()]);
            ind=dq.front();
        }
    }
    return;
}

vector<pair<unsigned int,unsigned int>> merge(vector<pair<unsigned int,unsigned int>> T){
    vector<pair<unsigned int,unsigned int>> ans;
    pair<unsigned int,unsigned int> temp=T[0];
    unsigned int i=1;
    unsigned int n=T.size();

    while(i<n){
        pair<unsigned int,unsigned int> curr=T[i];
        if (curr.first<=temp.second){
            if (curr.second>temp.second) temp.second=curr.second;
        }
        else{
            ans.push_back(temp);
            temp=curr;
        }
        i++;
    }
    ans.push_back(temp);
    return ans;
}


int main(){

    //ifstream wholeGenome("Genome_Anamika.txt");
    //string gen((istreambuf_iterator<char>(wholeGenome)),
    //                        istreambuf_iterator<char>());
    //gen.erase(remove(gen.begin(), gen.end(), '\n'), gen.end());    
    //wholeGenome.close(); 

    
    //int th=120;
    //long long int chunk=()/th;

    //vector<vector<long long int>> local_gen(th);
    //#pragma omp parallel num_threads(th)
    //{
    //    int thread_id=omp_get_thread_num();
    //    int start=thread_id*chunk;
    //    int end=(thread_id==th-1) ? gen.length():(thread_id+1)*chunk;
    //    
    //    if (start<gen.length()) local_gen[thread_id]=pre_hash(gen.substr(start,end-start));
    //}
    //vector<long long int> gen_hash;
    //for (int i=0; i<th; i++) gen_hash.insert(gen_hash.end(),local_gen[i].begin(),local_gen[i].end());
    
    //win_gen(gen_hash,25);

    
    // cout<<read<<endl;
    
    // omp_set_num_threads(50);
    

    // vector<vector<long long int>> local_read(th);
    // // gen_hash.reserve(gen.length());
    // #pragma omp parallel num_threads(th)
    // {
    //     int thread_id=omp_get_thread_num();
    //     int start=thread_id*chunk;
    //     int end=(thread_id==th-1) ? read.length():(thread_id+1)*chunk;
    //     // vector<long long int>local;
    //     //local.reserve(end-start);
    //     if (start<read.length()) local_read[thread_id]=pre_hash(read.substr(start,end-start));
        
    //         // #pragma omp critical
    //         // gen_hash.insert(gen_hash.end(),local.begin(),local.end());
        
    // }
    
    // ifstream inp1("g_len.txt"); 
    // long long int g_len=0;
    // inp1>>g_len;
    // inp1.close();
    
    ifstream inp2("g_hash.txt");
    unordered_map<unsigned int, vector<unsigned int>> g_hash;
    unsigned int g_len=0;
    inp2>>g_len;
    unsigned int key=0;
    while(inp2>>key){
    	auto vals=fileToVec<unsigned int>(inp2);
    	g_hash[key]=vals;
    }
    inp2.close();
    
    ifstream inp3("g_pos.txt");
    vector<unsigned int> g_pos(g_len);
    unsigned int idx=0, val_s=0;
    while(inp3>>val_s){
    	//inp3>>val_s;
        g_pos[idx++]=val_s;
    	//g_pos.push_back(make_pair(key_s,val_s));
    }
    inp3.close();
    
    cout<<g_len<<endl;
    cout<<g_hash.size()<<endl;
    cout<<g_pos.size()<<endl;
    // cout<<(sizeof(long long int) + sizeof(long long int) + sizeof(_Rb_tree_node_base)) * g_len + sizeof(_Rb_tree)<<endl;
    
    //window size=constant*read length
    //win_gen(gen_hash,25);
    
    ifstream readFile("read_ex.txt");
    string read((istreambuf_iterator<char>(readFile)),
                           istreambuf_iterator<char>());
    readFile.close();

    vector<unsigned int> read_hash=pre_hash(read);
    win_read(read_hash,25);
    
    unsigned int s=r.size();
    float e=0.1;
    float G=1/((2*(exp(e*16)))-1);

    cout<<"G="<<G<<endl;
    cout<<"s="<<s<<endl;
    cout<<read.length()<<endl;
    float m=s*G;
    m=ceil(m);
    cout<<"m="<<m<<endl;

    // cout<<"g_hash="<<sizeof(g_hash)<<endl;
    // cout<<"g_pos="<<sizeof(g_pos)<<endl;

    struct timeval start,end;
    gettimeofday(&start,NULL);
    ios_base::sync_with_stdio(false);

    vector<unsigned int> L;
    for (auto el: r){
        if (g_hash.find(el)!=g_hash.end()){
            for (auto pos:g_hash[el]) L.push_back(pos);
        }
    }
    sort(L.begin(),L.end());
    if (L.size()-m<0){
        cout<<L.size()-m<<endl;
        gettimeofday(&end,NULL);
        double time_taken;
        time_taken = (end.tv_sec - start.tv_sec) * 1e6;
        time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
        cout<<fixed<<time_taken<<setprecision(6)<<endl;
        return 0;
    }
    vector<pair<unsigned int,unsigned int>> T;

    unsigned int r_len=read.length();
    // cout<<r_len;
    for (unsigned int i=0; i<=L.size()-m; i++){
        unsigned int j=i+(m-1);
        if (L[j]-L[i]<r_len){
            if (L[j]-r_len+1<=0) T.push_back(make_pair(0,L[i]));
            else T.push_back(make_pair(L[j]-read.length()+1,L[i]));
        }
    } 
    if (T.size()==0){
        //cout<<L.size()-m<<endl;
        gettimeofday(&end,NULL);
        double time_taken;
        time_taken = (end.tv_sec - start.tv_sec) * 1e6;
        time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
        cout<<fixed<<time_taken<<setprecision(6)<<endl;
        return 0;
    }  

    // cout<<"T.size= "<<T.size()<<endl;
    // for (auto el:T) cout<<el.first<<"-"<<el.second<<endl; 

    vector<pair<unsigned int,unsigned int>> T_merged=merge(T);
    // cout<<"T_merged.size= "<<T_merged.size()<<endl;
    for (auto el: T_merged) cout<<el.first<<"-"<<el.second<<endl;
    
    vector<unsigned> ans;

    // vector<bool> check(T_merged[T_merged.size()-1].second+1,false);
    unordered_map<unsigned int,int> hash_map;
    for (auto el:r){
        hash_map[el]=0;
    } 

    unsigned int i=-1;
    unsigned int j=-1;

    //second condition

    sort(r.begin(),r.end());

    for (unsigned int k=0; k<T_merged.size(); k++){
        i=T_merged[k].first;
        j=T_merged[k].first+r_len-1;

        vector<int>read_map(r_len,0);

        int dirty=0;    
        unordered_map<unsigned int,int>hash_x_y=hash_map;
        //g_pos= pos->hash
        //hash from i to i+3.. then i+4 to i+7.
        
        unsigned int count=0;
        
        // for (int ind=i; ind<=j; ind++){
        //     long long int index=search_pair(g_pos,0,g_pos.size()-1,ind);
        //     if (index!=-1){
        //         long long int read_index=search(r,0,r.size()-1,g_pos[index].second);
        //         if (read_index!=-1){
        //             count++;
        //             read_map[ind-i]=1;
        //         }
        //     }
            
        // }
        //******************************


        // for (auto el:g_pos){
        //     if (el.first<=j && el.first>=i){
        //         if (hash_x_y.find(el.second)!=hash_x_y.end()){
        //             if (el.first==i) dirty=1;
        //             hash_x_y[el.second]=1;
        //             count++;
        //             }
        //         }
        //     }


        for (unsigned int ind=i; ind<=j; ind++){
            unsigned int h=g_pos[ind];
            if (hash_x_y.find(h)!=hash_x_y.end()){
                if (ind==i) dirty=1;
                hash_x_y[h]=1;
                count++;
            }
                            
        }   


        float jac=float(count)/float(s);
        if (jac>=G){
        //cout<<"index "<<i<<"->"<<jac<<endl;
            ans.push_back(i);
            }

        while(i<=T_merged[k].second && j<=g_len-2){    

            // if (read_map[0]==1) count--;

            // //left shift
            // int carry=0;
            // for (int ind=0; ind<read_map.size()-1; ++ind) {
            //         int temp = read_map[i];
            //         read_map[i] = (read_map[i] << 1) | carry;
            //         carry = temp >> (sizeof(int) * 8 - 1);
            //     }
            
            // //for j+1
            // int index=search_pair(g_pos,0,g_pos.size()-1,j+1);
            // if (index!=-1){
            //     int read_index=search(r,0,r.size()-1,g_pos[index].second);
            //     if (read_index!=-1){
            //         count++;
            //         read_map[read_map.size()-1]=1;
            //     }
            //     else read_map[read_map.size()-1]=0;
            // }
            // else read_map[read_map.size()-1]=0;


            if (dirty==1) count--;
            
            unsigned int h1=g_pos[i+1];
            if (hash_x_y.find(h1)!=hash_x_y.end()){
                if (hash_x_y[h1]==1) dirty=1;
                }
                else{
                    dirty=0;
                }
                 


            unsigned int h2=g_pos[j+1];
            if (hash_x_y.find(h2)!=hash_x_y.end()){
                if (hash_x_y[h2]==0){
                    count++;
                    hash_x_y[h2]=1;
                }
            }
                    
            



            // for (auto el:g_pos){
            //     if (el.first==i+1){
            //         if (hash_x_y.find(el.second)!=hash_x_y.end()){
            //             if (hash_x_y[el.second]==1) dirty=1;
            //         }
            //         else{
            //             dirty=0;
            //         }
            //     }

            //     else if (el.first==j+1){
            //         if (hash_x_y.find(el.second)!=hash_x_y.end()){
            //             if (hash_x_y[el.second]==0){
            //                 count++;
            //                 hash_x_y[el.second]=1;
            //             }
            //         }
            //     }
            // }

            float jac=float(count)/float(s);
            if (jac>=G){
            //cout<<"index "<<i+1<<"->"<<jac<<endl;
            ans.push_back(i+1);
            }
                i++; j++;
        }
        cout<<endl;
    }
    
    cout<<ans.size()<<endl;
    
    gettimeofday(&end,NULL);
    double time_taken;
    time_taken = (end.tv_sec - start.tv_sec) * 1e6;
    time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
    cout<<fixed<<time_taken<<setprecision(6)<<endl;
    return 0;
    

    // for (auto el: ans) cout<<el<<" ";
    // cout<<endl;

// try with a read of huge length(5000...) picked up from the genome itself
    // auto stop = high_resolution_clock::now();
    // auto duration = duration_cast<seconds>(stop - start);
    // cout << duration.count() << endl;
    
}

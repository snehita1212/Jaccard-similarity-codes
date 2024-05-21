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

vector<unsigned int> win_read(vector<unsigned int> hash, unsigned int k){
    vector<unsigned int> r;
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
    return r;
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

vector<unsigned int> mapping_one(string read, float G, unordered_map<unsigned int, vector<unsigned int>> g_hash, vector<unsigned int> g_pos, unsigned int g_len){
    vector<unsigned int> read_hash=pre_hash(read);
    vector<unsigned int> r=win_read(read_hash,25);
    vector<unsigned> ans;

    unsigned int s=r.size();
    //cout<<"s="<<s<<endl;
    //cout<<read.length()<<endl;
    float m=s*G;
    m=ceil(m);
    //cout<<"m="<<m<<endl;

    //pre-filter
    vector<unsigned int> L;
    for (auto el: r){
        if (g_hash.find(el)!=g_hash.end()){
            for (auto pos:g_hash[el]) L.push_back(pos);
        }
    }
    sort(L.begin(),L.end());
    if (L.size()-m<0){
        return ans;
    }
    vector<pair<unsigned int,unsigned int>> T;
    unsigned int r_len=read.length();
    for (unsigned int i=0; i<=L.size()-m; i++){
        unsigned int j=i+(m-1);
        if (L[j]-L[i]<r_len){
            if (L[j]-r_len+1<=0) T.push_back(make_pair(0,L[i]));
            else T.push_back(make_pair(L[j]-read.length()+1,L[i]));
        }
    } 
    if (T.size()==0){
        return ans;
    }
    // cout<<"T.size= "<<T.size()<<endl;
    // for (auto el:T) cout<<el.first<<"-"<<el.second<<endl;
    vector<pair<unsigned int,unsigned int>> T_merged=merge(T);
    // cout<<"T_merged.size= "<<T_merged.size()<<endl;
    //for (auto el: T_merged) cout<<el.first<<"-"<<el.second<<endl;
    
    //second condition
    
    unordered_map<unsigned int,int> hash_map;
    for (auto el:r){
        hash_map[el]=0;
    }
    unsigned int i=-1;
    unsigned int j=-1;
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
            continue;
            }

        while(i<=T_merged[k].second && j<=g_len-2){    
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

            float jac=float(count)/float(s);
            if (jac>=G){
            //cout<<"index "<<i+1<<"->"<<jac<<endl;
            ans.push_back(i+1);
            break;
            }
                i++; j++;
        }
        // cout<<endl;
    }
    // cout<<ans.size()<<endl;
    return ans;   
}

void writeSAM(const string& filename, vector<pair<string, vector<unsigned int>>> alignmentResults) {
    ofstream samFile(filename);

    // Write SAM header
    samFile << "@HD\tVN:1.6\tSO:coordinate\n";
    samFile << "@SQ\tSN:U00096.3\tLN:4641652\n"; // Example reference sequence, adjust accordingly

    // Write alignment results to SAM file
    for (auto result : alignmentResults) {
        string fullreadID = result.first;
        istringstream iss(fullreadID);
        string readID;
        iss>>readID;
        long long int len=result.second.back();
        result.second.pop_back();
        for (auto ind: result.second){ 
            long long int pos=ind;
        // const std::string& refName = std::get<1>(result);
        // int pos = std::get<2>(result);

        // Write the SAM alignment line
        samFile << readID << "\t"   // QNAME
                << "0" << "\t"      // FLAG (assuming single-end read with no specific flags)
                << len << "\t"
                << "U00096.3" << "\t"  // RNAME
                << pos << "\t"      // POS
                << "255" << "\t"    // MAPQ (dummy value)
                << "NA" << "\t"    // CIGAR (assuming full match for simplicity)
                << "*" << "\t"      // RNEXT
                << "0" << "\t"      // PNEXT
                << "0" << "\t"      // TLEN
                << "*" << "\t"      // SEQ (unknown sequence)
                << "*" << "\n";     // QUAL (unknown quality)
        }   
    }

    samFile.close();
}

int main(){
    
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
        g_pos[idx++]=val_s;
        //g_pos.push_back(make_pair(key_s,val_s));
    }
    inp3.close();
    
    // cout<<g_len<<endl;
    // cout<<g_hash.size()<<endl;
    // cout<<g_pos.size()<<endl;
    // cout<<(sizeof(long long int) + sizeof(long long int) + sizeof(_Rb_tree_node_base)) * g_len + sizeof(_Rb_tree)<<endl;
    
    //win_gen(gen_hash,25);
   
    ifstream readFile("ERR9979375.fastq");
    string read;

    unsigned int step=1;
    vector<pair<string,string>> reads;
    string id, r;
    while(getline(readFile,read)){
        if (step%4==1){
            id=read;
            step++;
        }
        else if (step%4==2){
            r=read;
            reads.push_back(make_pair(id,r));
            step++;
        }
        else if (step%4==3){
            step++;
            continue;
        }
        else{
            step=1;
            // step++;
            continue;
        }
    }
    readFile.close();
    
    float e=0.1;
    float G=1/((2*(exp(e*16)))-1);

    struct timeval start,end;
    gettimeofday(&start,NULL);
    ios_base::sync_with_stdio(false);
    

    // int n_threads=20000;
    vector<pair<string, vector<unsigned int>>> res(reads.size());
    // unsigned int chunk=reads.size()/n_threads;
    
    // omp_set_num_threads(10000);
    #pragma omp parallel for
    for (size_t i = 0; i < reads.size(); ++i){
        auto& item = reads[i];
        vector<unsigned int> result=mapping_one(item.second,G,g_hash,g_pos,g_len);
        result.push_back(item.second.length());
        
        res[i]=make_pair(reads[i].first,result);
        }      
    

    
    // for (int i=0; i<n_threads; i++) res.insert(res.end(),local_res[i].begin(),local_res[i].end());

    gettimeofday(&end,NULL);
    double time_taken;
    time_taken = (end.tv_sec - start.tv_sec) * 1e6;
    time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
    cout<<fixed<<time_taken<<setprecision(6)<<endl;

    string refName="U00096.3";

    writeSAM("output.sam",res);


    
    return 0;

    //second condition   
}


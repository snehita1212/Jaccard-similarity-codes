#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/time.h>
#include <omp.h>
#include <pthread.h>
#include <sys/ipc.h>
#include <sys/shm.h>
//# include "shared.h"

// #include <chrono>
using namespace std;

unordered_map<unsigned int,vector<unsigned int>> g_hash;
vector<unsigned int> g_pos; 

template<typename T>
void vecToFile(ofstream& out, vector<T>& vec){
	out<<vec.size()<<" ";
	for (auto& ele: vec) out<<ele<<" ";
	out<<endl;
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
    for (unsigned int i=0; i<kgram.size(); i++){
        unsigned int num=0;
        for (char c:kgram[i]){
            num=num*2+(c-'0');
        }
        
        hash.push_back(num);
    }
    
    return hash;
}

void win_gen(vector<unsigned int> hash, unsigned int k){
    deque<unsigned int> dq;
    
    for (unsigned int i=0; i<k; i++){
        while(!dq.empty() && hash[i]<=hash[dq.back()]) dq.pop_back();
        g_pos.push_back(hash[i]);
        
        dq.push_back(i);
    }

    g_hash[hash[dq.back()]].push_back(dq.back());
    // res[hash[dq.back()]].push_back(dq.back());

    //g_pos.push_back(make_pair(dq.back(),hash[dq.back()]));
    // g_pos[dq.back()]=hash[dq.back()];

    unsigned int ind=dq.front();
    for (unsigned int i=k; i<hash.size(); i++){
        g_pos.push_back(hash[i]);
        
        if (dq.front()==i-k) dq.pop_front();
        while(!dq.empty() && hash[i]<=hash[dq.back()]) dq.pop_back();
        dq.push_back(i);
        if (ind!=dq[0]){
            g_hash[hash[dq.front()]].push_back(dq.front());
            // res.push_back({hash[dq.front()],dq.front()});
            
            // g_pos[dq.front()]=hash[dq.front()];
            //g_pos.push_back(make_pair(dq.front(),hash[dq.front()]));

            ind=dq.front();
        }
    }
    return;
}


int main(){

    ifstream wholeGenome("Genome_Anamika.txt");
    string gen((istreambuf_iterator<char>(wholeGenome)),
                            istreambuf_iterator<char>());
    gen.erase(remove(gen.begin(), gen.end(), '\n'), gen.end());    
    wholeGenome.close(); 
    
    struct timeval start,end;
    gettimeofday(&start,NULL);
    ios_base::sync_with_stdio(false);
    
    int th=50;
    unsigned int chunk=gen.length()/th;

    vector<vector<unsigned int>> local_gen(th);
    
    #pragma omp parallel num_threads(th)
    {
        int thread_id=omp_get_thread_num();
       int start=thread_id*chunk;
        int end=(thread_id==th-1) ? gen.length():(thread_id+1)*chunk;
        
        if (start<gen.length()) local_gen[thread_id]=pre_hash(gen.substr(start,end-start));
       
        
    }
    vector<unsigned int> gen_hash;
    for (int i=0; i<th; i++) gen_hash.insert(gen_hash.end(),local_gen[i].begin(),local_gen[i].end());
    
    win_gen(gen_hash,25);
    //sort(g_pos.begin(),g_pos.end());

    gettimeofday(&end,NULL);
    double time_taken;
    time_taken = (end.tv_sec - start.tv_sec) * 1e6;
    time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
    cout<<fixed<<time_taken<<setprecision(6)<<endl;
    
    cout<<gen_hash.size()<<endl;
    cout<<g_hash.size()<<endl;
    cout<<g_pos.size()<<endl;
    //unordered_map<long long int,vector<long long int>> g_hash=win_gen_hash(gen_hash,25);
    //unordered_map<long long int, long long int> g_pos=win_gen_pos(gen_hash,25);
    
    unsigned int g_len=gen.length();

    //read to genome_output.txt
    
    // ofstream out1("gen_len.txt");
    // 	out1<<g_len<<endl;
    // out1.close();
    
    	//output<<"#\n";
    ofstream out2("g_hash.txt");
    	//g_hash(int->vector)
        out2<<g_len<<endl;
    	for (auto pair: g_hash){
    		out2<<pair.first<<" ";
    		vecToFile(out2,pair.second);
    	}
    out2.close();
    
    ofstream out3("g_pos.txt");
    	//g_pos(int->int)
    	for (int i=0; i<g_pos.size(); i++){
    		out3<<g_pos[i]<<endl;
    	}
    out3.close();
    
    //else return 1;
    
    
    return 0;

    //return 0;
}
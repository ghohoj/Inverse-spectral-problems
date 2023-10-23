/*
 * @Description: 请加入备注来说明这个文件的意义
 * @Author: catgod
 * @Date: 2023-10-20 18:09:39
 * @LastEditTime: 2023-10-21 21:05:02
 * @FilePath: /Inverse spectral problems/test.cpp
 */
/*
 *                        _oo0oo_
 *                       o8888888o
 *                       88" . "88
 *                       (| -_- |)
 *                       0\  =  /0
 *                     ___/`---'\___
 *                   .' \\|     |// '.
 *                  / \\|||  :  |||// \
 *                 / _||||| -:- |||||- \
 *                |   | \\\  - /// |   |
 *                | \_|  ''\---/''  |_/ |
 *                \  .-\__  '-'  ___/-. /
 *              ___'. .'  /--.--\  `. .'___
 *           ."" '<  `.___\_<|>_/___.' >' "".
 *          | | :  `- \`.;`\ _ /`;.`/ - ` : | |
 *          \  \ `_.   \_ __\ /__ _/   .-` /  /
 *      =====`-.____`.___ \_____/___.-`___.-'=====
 *                        `=---='
 * 
 * 
 *      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 *            佛祖保佑       永不宕机     永无BUG
 */

#include"MatrixNumerovmethod/def.h"
#include"MatrixNumerovmethod/Numerov.h"
#include"MatrixNumerovmethod/rV.h"
#include"MatrixNumerovmethod/ReadandWrite.h"



void test1(){
    vector<vector<Real>> qs;
    vector<vector<Real>> lams;
    for(int j=0;j<10;j++){
        vector<Real> v=ployV();
        vector<Real> eigval;
        MatrixNumerov(v,eigval,20);
        qs.push_back(v);
        lams.push_back(eigval);

        cout<<endl;
        cout<<"第"<<j<<"轮结果:";
        for(auto i:eigval){ 
            cout<<i<<" ";
        }
        cout<<endl;
    }
    write(qs,lams);
}

void test2(){
    int n=4;
    string address="./data/train.json";
    fstream f;
    f.open(address,ios::out);
    f<<"["<<endl;
    int trainnum=1000;
    for(int j=0;j<trainnum;j++){
        cout<<"第"<<j<<"轮 ";
        vector<Real> lam;

        vector<Real> cof(4);
        for(int i=0;i<n;i++){
            cof[i]=(myrand()-0.5)*10;
        }
        vector<Real> tmp=CosV(cof);
        MatrixNumerov(tmp,lam,10);

        f<<"{";
        f<<" \"cof\":[ ";
        for(int i=0;i<n-1;i++){
            f<<setprecision(10)<<cof[i]<<",";
        }
        f<<setprecision(10)<<cof[n-1];
        f<<"],";
        f<<" \"lambda\":[ ";
        sort(lam.begin(),lam.end());
        for(int i=0;i<lam.size()-1;i++){
            f<<setprecision(10)<<lam[i]<<",";
        }
        f<<setprecision(10)<<lam[lam.size()-1];
        f<<"]";
        f<<"}";
        if(j!=trainnum-1){
            f<<","<<endl;
        }
        else{
            f<<endl;
        }
        cout<<endl;
    }
    f<<"]"<<endl;
    f.close();
}


void test3(){
    for(int k=0;k<10;k++){
        int n=4;
        string address="./data/test"+to_string(k)+".json";
        fstream f;
        f.open(address,ios::out);
        f<<"["<<endl;
        int trainnum=100;
        for(int j=0;j<trainnum;j++){
            cout<<"第"<<j<<"轮 ";
            vector<Real> lam;
            vector<Real> cof(4) ;
            for(int i=0;i<n;i++){
                cof[i]=(myrand()-0.5)*10;
            }
            vector<Real> tmp=CosV(cof);
            MatrixNumerov(tmp,lam,10);

            f<<"{";
            f<<" \"cof\":[ ";
            for(int i=0;i<n-1;i++){
                f<<setprecision(10)<<cof[i]<<",";
            }
            f<<setprecision(10)<<cof[n-1];
            f<<"],";
            f<<" \"lambda\":[ ";
            sort(lam.begin(),lam.end());
            for(int i=0;i<lam.size()-1;i++){
                f<<setprecision(10)<<lam[i]<<",";
            }
            f<<setprecision(10)<<lam[lam.size()-1];
            f<<"]";
            f<<"}";
            if(j!=trainnum-1){
                f<<","<<endl;
            }
            else{
                f<<endl;
            }
            cout<<endl;
        }
        f<<"]"<<endl;
        f.close();
    }
    
} 


int main()
{
    test3();
}
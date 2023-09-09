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
    auto v=testV();
    vector<Real> eigval;
    MatrixNumerov(v,eigval,4);
    for(auto i:eigval){
        cout<<i<<endl;
    }
}

void test2(){
    vec q(precise);
    vector<Real> lam;

    for(int i=0;i<precise;i++){
        q.insert(i)=0;
    }
    for(int i=0;i<20;i++){
        lam.push_back((i+1)*(i+1)*M_PI*M_PI);
    }
    write(q,lam);
} 

int main()
{
    test2();
}
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include<vector>
#include<stack>
#include <cmath>

using namespace std;

double pi = 2*acos(0.0);

double eyeX,eyeY,eyeZ,lookX,lookY,lookZ,upX,upY,upZ,fovY,aspectRatio,near,far;
vector<vector<double>> transformationMatrix( 4 , vector<double> (4, 0));
stack<vector<vector<double>>> transformationStates;
vector<vector<vector<double>>> stage1Triangles,stage2Triangles;

ifstream sceneFile;
ofstream stage1File;
ofstream stage2File;
ofstream stage3File;

void initializeTransformationMatrix(){
    for(int i=0;i<4;i++)
        transformationMatrix[i][i]=1;
    transformationStates.push(transformationMatrix);
}

void viewTriangle(vector<vector<double>> triangleMatrix){
    for(int i=0;i<4;i++){
        for(int j=0;j<3;j++)
            cout<<fixed<<setprecision(7)<<triangleMatrix[i][j]<<"\t";
        cout<<endl;
    }
}

void viewTransformation(vector<vector<double>> transformationMatrix){
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++)
            cout<<fixed<<setprecision(7)<<transformationMatrix[i][j]<<"\t";
        cout<<endl;
    }
}

vector<vector<double>> transform(vector<vector<double>>transformMatrix,vector<vector<double>> operandMatrix){
    vector<vector<double>> outputMatrix( transformMatrix.size() , vector<double> (operandMatrix[0].size(), 1));
    if(transformMatrix[0].size()==operandMatrix.size())
        for(int i=0;i<transformMatrix.size();i++)
            for(int j=0;j<operandMatrix[0].size();j++){
                double sum=0;
                for(int k=0;k<operandMatrix.size();k++)
                    sum+=transformMatrix[i][k]*operandMatrix[k][j];
                outputMatrix[i][j]=sum;
            }
    return outputMatrix;
}

vector<vector<double>> copyTransformation(vector<vector<double>>state){
    vector<vector<double>> outputMatrix( 4 , vector<double> (4, 1));
    for(int i=0;i<4;i++)
        for(int j=0;j<4;j++)
            outputMatrix[i][j]=state[i][j];
    return outputMatrix;
}

vector<double> vectorScalarMultiple1D(vector<double>x, double m){
    vector<double> outputVector(x.size(), 0);
    for(int i=0;i<x.size();i++)
        outputVector[i]=m*x[i];
    return outputVector;
}

vector<double> getVector(double x,double y,double z){
    vector<double> outputVector(3, 0);
    outputVector[0]=x;
    outputVector[1]=y;
    outputVector[2]=z;
    return outputVector;
}

double vectorDotMultiple1D(vector<double>v1, vector<double>v2){
    double sum=0;
    for(int i=0;i<v1.size();i++)
        sum+=v1[i]*v2[i];
    return sum;
}

vector<double> vectorAddition1D(vector<double>v1, vector<double>v2){
    vector<double> outputVector(v1.size(), 0);
    for(int i=0;i<v1.size();i++)
        outputVector[i]=v1[i]+v2[i];
    return outputVector;
}

vector<double> vectorCrossMultiple1D(vector<double>v1, vector<double>v2){
    vector<double> outputVector(v1.size(), 0);
    outputVector[0]=v1[1]*v2[2]-v1[2]*v2[1];
    outputVector[1]=v1[2]*v2[0]-v1[0]*v2[2];
    outputVector[2]=v1[0]*v2[1]-v1[1]*v2[0];
    return outputVector;
}

vector<double> getNormalizedVector1D(vector<double>v){
    double squareSum=0;
    for(int i=0;i<v.size();i++)
        squareSum+=pow(v[i],2);
    double value=sqrt(squareSum);
    return vectorScalarMultiple1D(v,1/value);
}

vector<double> getDirectionalVector(int direction){
    vector<double> outputVector(3, 0);
    outputVector[direction]=1;
    return outputVector;
}

//𝑅(𝑥⃗, 𝑎⃗, 𝜃) = cos 𝜃 𝑥⃗ + (1 − cos 𝜃)(𝑎⃗ ∙ 𝑥⃗)𝑎⃗ + sin 𝜃 (𝑎⃗ × 𝑥⃗)
vector<double> rodriguesFunction(vector<double> x,vector<double> a, double thetaRadian){
    vector<double> t_1=vectorScalarMultiple1D(x,cos(thetaRadian));
    vector<double> t_2=vectorScalarMultiple1D(a,(1-cos(thetaRadian))*vectorDotMultiple1D(a,x));
    vector<double> t_3=vectorScalarMultiple1D(vectorCrossMultiple1D(a,x),sin(thetaRadian));
    return vectorAddition1D(t_1,vectorAddition1D(t_2,t_3));
}

void writeTriangle(ofstream &file,vector<vector<double>> triangle){
        for(int i=0;i<3;i++)
            file<<fixed<<setprecision(7)<<triangle[0][i]<<"\t"<<triangle[1][i]<<"\t"<<triangle[2][i]<<endl;
        file<<endl;
}

void takeInput(){
    
    if (sceneFile.is_open()){ 
        sceneFile>>eyeX>>eyeY>>eyeZ>>lookX>>lookY>>lookZ>>upX>>upY>>upZ>>fovY>>aspectRatio>>near>>far;
        string input;
        while(true){
            sceneFile>>input;
            if(input=="end")break;
            if(input=="triangle"){
                vector<vector<double>> triangleMatrix( 4 , vector<double> (3, 1));
                for(int i=0;i<3;i++)
                    sceneFile>>triangleMatrix[0][i]>>triangleMatrix[1][i]>>triangleMatrix[2][i];
                vector<vector<double>> transformedTriangleMatrix=transform(transformationStates.top(),triangleMatrix);
                stage1Triangles.push_back(transformedTriangleMatrix);
                writeTriangle(stage1File,transformedTriangleMatrix);
            }else if(input=="translate"){
                vector<vector<double>> translateMatrix( 4 , vector<double> (4, 0));
                for(int i=0;i<4;i++)
                    translateMatrix[i][i]=1;
                for(int i=0;i<3;i++)
                    sceneFile>>translateMatrix[i][3];
                vector<vector<double>> newTransformedMatrix=transform(transformationStates.top(),translateMatrix);
                transformationStates.pop();
                transformationStates.push(newTransformedMatrix);
            }
            else if(input=="scale"){
                vector<vector<double>> scaleMatrix( 4 , vector<double> (4, 0));
                for(int i=0;i<4;i++)
                    scaleMatrix[i][i]=1;
                for(int i=0;i<3;i++)
                    sceneFile>>scaleMatrix[i][i];
                vector<vector<double>> newTransformedMatrix=transform(transformationStates.top(),scaleMatrix);
                transformationStates.pop();
                transformationStates.push(newTransformedMatrix);
            }
            else if(input=="rotate"){
                vector<double> rotationVector( 3 , 0);
                double angleDegree;
                sceneFile>>angleDegree;
                double angleRadian=angleDegree*pi/180.0;
                for(int i=0;i<3;i++)
                    sceneFile>>rotationVector[i];
                
                rotationVector=getNormalizedVector1D(rotationVector);
                vector<double> c1=rodriguesFunction(getDirectionalVector(0),rotationVector,angleRadian);
                vector<double> c2=rodriguesFunction(getDirectionalVector(1),rotationVector,angleRadian);
                vector<double> c3=rodriguesFunction(getDirectionalVector(2),rotationVector,angleRadian);

                vector<vector<double>> rotationMatrix( 4 , vector<double> (4, 0));
                for(int i=0;i<4;i++)
                    rotationMatrix[i][i]=1;

                for(int i=0;i<3;i++){
                    rotationMatrix[i][0]=c1[i];
                    rotationMatrix[i][1]=c2[i];
                    rotationMatrix[i][2]=c3[i];
                }

                vector<vector<double>> newTransformedMatrix=transform(transformationStates.top(),rotationMatrix);
                transformationStates.pop();
                transformationStates.push(newTransformedMatrix);
            }
            else if(input=="push"){
                vector<vector<double>> newTransformedMatrix=copyTransformation(transformationStates.top());
                transformationStates.push(newTransformedMatrix);
            }
            else if(input=="pop")
                transformationStates.pop();
            
        }
        sceneFile.close(); 
        stage1File.close(); 
    }
}

void executeStage2(){
    vector<double> l=vectorAddition1D(getVector(lookX,lookY,lookZ),vectorScalarMultiple1D(getVector(eyeX,eyeY,eyeZ),-1));
    l=getNormalizedVector1D(l);
    vector<double> r=vectorCrossMultiple1D(l,getVector(upX,upY,upZ));
    r=getNormalizedVector1D(r);
    vector<double> u=vectorCrossMultiple1D(r,l);

    vector<vector<double>> tMatrix( 4 , vector<double> (4, 0));
    for(int i=0;i<4;i++)
        tMatrix[i][i]=1;

    tMatrix[0][3]=-eyeX;
    tMatrix[1][3]=-eyeY;
    tMatrix[2][3]=-eyeZ;

    vector<vector<double>> rMatrix( 4 , vector<double> (4, 0));
    for(int i=0;i<4;i++)
        tMatrix[i][i]=1;
    
    for(int i=0;i<3;i++){
        rMatrix[0][i]=r[i];
        rMatrix[1][i]=u[i];
        rMatrix[2][i]=-l[i];
    }

    vector<vector<double>> vMatrix=transform(rMatrix,tMatrix);

    for(int i=0;i<stage1Triangles.size();i++){
        vector<vector<double>> transformedTriangleMatrix=transform(vMatrix,stage1Triangles[i]);
        stage2Triangles.push_back(transformedTriangleMatrix);
        writeTriangle(stage2File,transformedTriangleMatrix);
    }

    stage2File.close(); 
}


void executeStage3(){
    double fovX=fovY*aspectRatio;
    double t=near*tan(fovY/2);
    double r=near*tan(fovX/2);

    vector<vector<double>> pMatrix( 4 , vector<double> (4, 0));
    pMatrix[0][0]=near/r;
    pMatrix[1][1]=near/t;
    pMatrix[2][2]=-(far+near)/(far-near);
    pMatrix[2][3]=-2*far*near/(far-near);
    pMatrix[3][2]=-1;

    // viewTransformation(pMatrix);

    for(int i=0;i<stage2Triangles.size();i++){
        vector<vector<double>> transformedTriangleMatrix=transform(pMatrix,stage2Triangles[i]);
        //stage2Triangles.push_back(transformedTriangleMatrix);
        writeTriangle(stage3File,transformedTriangleMatrix);
    }

    stage3File.close(); 


}

int main(){

    initializeTransformationMatrix();

    sceneFile.open("Test Cases/1/scene.txt",ios::in); 
    stage1File.open("outputs/1/stage1.txt"); 
    stage2File.open("outputs/1/stage2.txt");
    stage3File.open("outputs/1/stage3.txt");

    //stage1
    takeInput();    

    //stage2
    executeStage2();

    //stage3
    executeStage3();

    return 0;
}
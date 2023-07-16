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

vector<vector<double>> transformTriangle(vector<vector<double>>transformMatrix,vector<vector<double>> triangleMatrix){
    vector<vector<double>> outputMatrix( 4 , vector<double> (3, 1));
    for(int i=0;i<4;i++)
        for(int j=0;j<3;j++){
            double sum=0;
            for(int k=0;k<4;k++)
                sum+=transformMatrix[i][k]*triangleMatrix[k][j];
            outputMatrix[i][j]=sum;
        }
    return outputMatrix;
}

vector<vector<double>> transformState(vector<vector<double>>state1,vector<vector<double>> state2){
    vector<vector<double>> outputMatrix( 4 , vector<double> (4, 1));
    for(int i=0;i<4;i++)
        for(int j=0;j<4;j++){
            double sum=0;
            for(int k=0;k<4;k++)
                sum+=state1[i][k]*state2[k][j];
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

void takeInput(){
    ifstream sceneFile;
    ofstream stage1File;
    sceneFile.open("Test Cases/4/scene.txt",ios::in); 
    stage1File.open("outputs/4/stage1.txt"); 


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
                vector<vector<double>> transformedTriangleMatrix=transformTriangle(transformationStates.top(),triangleMatrix);
                //viewTriangle(transformedTriangleMatrix);
                for(int i=0;i<3;i++)
                    stage1File<<fixed<<setprecision(7)<<transformedTriangleMatrix[0][i]<<"\t"<<transformedTriangleMatrix[1][i]<<"\t"<<transformedTriangleMatrix[2][i]<<endl;
                stage1File<<endl;
            }else if(input=="translate"){
                vector<vector<double>> translateMatrix( 4 , vector<double> (4, 0));
                for(int i=0;i<4;i++)
                    translateMatrix[i][i]=1;
                for(int i=0;i<3;i++)
                    sceneFile>>translateMatrix[i][3];
                vector<vector<double>> newTransformedMatrix=transformState(transformationStates.top(),translateMatrix);
                transformationStates.pop();
                transformationStates.push(newTransformedMatrix);
            }
            else if(input=="scale"){
                vector<vector<double>> scaleMatrix( 4 , vector<double> (4, 0));
                for(int i=0;i<4;i++)
                    scaleMatrix[i][i]=1;
                for(int i=0;i<3;i++)
                    sceneFile>>scaleMatrix[i][i];
                vector<vector<double>> newTransformedMatrix=transformState(transformationStates.top(),scaleMatrix);
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

                vector<vector<double>> newTransformedMatrix=transformState(transformationStates.top(),rotationMatrix);
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

int main(){

    initializeTransformationMatrix();
    takeInput();    

    return 0;
}
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 19:22:24 2016

@author: dinos
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import networkx as nx



R=5;
space=1
xyDim = 5
mobile_nodes=2
time_instances=4
positions=np.zeros(((time_instances,mobile_nodes,2)))
var= np.zeros((mobile_nodes))
VarLimit=1
F=np.identity(2)
Delta=1


def AnchorNodes( space,xyDim ):
    '''
    space represents the number of spaces in between the number of base nodes
    on the x or the y axis
    '''
    b=np.arange(0,space+1,1)*xyDim/space
    X,Y=np.meshgrid(b,b)
    X = np.reshape(X,(len(b)**2,1))
    Y= np.reshape(Y,(len(b)**2,1))
    m_extende = np.hstack([X, Y])
    '''
    the distance in between the nodes which are one number less than those on 
    the x axis is half spaced compared to those on the x axis
    '''
    c=np.arange(1,space*2,2)*xyDim/space/2
    X1,Y1=np.meshgrid(c,c)
    X1 = np.reshape(X1,(len(c)**2,1))
    Y1= np.reshape(Y1,(len(c)**2,1))
    m_extend = np.hstack([X1, Y1])
    anchors=np.vstack((m_extende,m_extend))
    return anchors
    
anchors=(AnchorNodes(space,xyDim))




for i in range(mobile_nodes):
    positions[0,i,:]=[random.uniform(0, xyDim),random.uniform(0, xyDim)]
    var[i]=random.uniform(1, VarLimit)
    


for k in range(time_instances-1):
        '''
        plt.figure()
        plt.scatter(anchors[:,0], anchors[:,1],marker='s',color=['green'])
        plt.scatter(positions[k,:,0], positions[k,:,1],color=['red'])
        plt.show()
        '''
        for i in range(mobile_nodes):
            positions[k+1,i,:]=np.dot(positions[k,i,:],F)+np.dot(Delta*var[i],np.dot(np.random.rand(1,2),np.identity(2)))
            if positions[k+1,i,0]>(xyDim):
                positions[k+1,i,:] =[positions[k+1,i,0]-xyDim,positions[k+1,i,1]];
            if positions[k+1,i,0]<0:
                positions[k+1,i,:] =[positions[k+1,i,0]+xyDim,positions[k+1,i,1]];
            if positions[k+1,i,1]>(xyDim):
                positions[k+1,i,:] =[positions[k+1,i,0],positions[k+1,i,1]-xyDim];
            if positions[k+1,i,1]<0:
                positions[k+1,i,:] =[positions[k+1,i,0],positions[k+1,i,1]+xyDim];
'''       
plt.figure()
plt.scatter(anchors[:,0], anchors[:,1],marker='s',color=['green'])
plt.scatter(positions[k+1,:,0], positions[k+1,:,1],color=['red'])
plt.show()
'''


Nodes=np.zeros(((time_instances,mobile_nodes+len(anchors),2)));
for i in range(time_instances):
    Nodes[i,:,:]= np.vstack(( positions[i,:,:],anchors));   


Distance=np.zeros(((time_instances,mobile_nodes,mobile_nodes+len(anchors))));
Probabilities=np.zeros(((time_instances,mobile_nodes,mobile_nodes+len(anchors))));
Connectivity_status=np.zeros(((time_instances,len(anchors)+mobile_nodes,len(anchors)+mobile_nodes)));


alpha=0.6;
for L in range(time_instances):
    k=0;
    for i in range(mobile_nodes):
        k=k+1;
        for j in range(k,mobile_nodes+len(anchors)):
            Distance[L,i,j]=np.linalg.norm(positions[L,i,:]-Nodes[L,j,:])+np.random.rand(1);        
            Probabilities[L,i,j]=np.exp(-Distance[L,i,j]/(2*R^2));


'''
Initialize the Markov chain 
'''
Connectivity_status[0,0:mobile_nodes,:]=Probabilities[0,0:mobile_nodes,:]> np.random.uniform(0,1,(mobile_nodes,len(Probabilities[0,0])));           
for L in range(1,time_instances):
    k=0;
    for i in range(mobile_nodes):
        k=k+1;
        for j in range(k,mobile_nodes+len(anchors)):
            beta=(Probabilities[L,i,j]*alpha)/(1-Probabilities[L,i,j]);
            K=np.hstack([1-alpha, alpha]);
            M=np.hstack([beta, 1-beta]);
            Transition=np.vstack((K,M));
            if Connectivity_status[L-1,i,j]==1:
                Connectivity_status[L,i,j]=Transition[0,0]>np.random.uniform(0,1)
            else:
                Connectivity_status[L,i,j]=Transition[0,1]>np.random.uniform(0,1)  
    colors=[]; 
    for n in range(len(anchors)+mobile_nodes):
        if n < mobile_nodes:
            colors.append('g')
        else:
            colors.append('b')       
    Gr=nx.from_numpy_matrix( Connectivity_status[L,:,:]);
    pos=Nodes[L,:,:];
    plt.figure()
    nx.draw(Gr,pos,with_labels=True,node_color=colors) 
    '''
    plt.scatter(anchors[:,0], anchors[:,1],marker='s',s=500,color=['green'])
    plt.scatter(positions[L,:,0], positions[L,:,1],s=500,color=['red'])
    '''
    plt.show()
    
function [Distance, Probabilities,Anchors,Nodes,positions] = Data(R,space,Var, time_instances,mobile_nodes,xyDim)
b=(0:space)*xyDim./space;
[X,Y]=meshgrid(b,b);
X=reshape(X,length(b)^2,1);
Y=reshape(Y,length(b)^2,1);
First=[X,Y];


c=(1:2:2*space-1)*xyDim./(2*space);
[X1,Y1]=meshgrid(c,c);
X1=reshape(X1,length(c)^2,1);
Y1=reshape(Y1,length(c)^2,1);
Second=[X1,Y1];
Anchors=[First;Second];

F=eye(2);
Delta=1;
var=zeros(1,mobile_nodes);
positions=zeros(mobile_nodes,2,time_instances+1);
for i=1:mobile_nodes
    positions(i,:,1)=[randi([0 xyDim],1),randi([0 xyDim],1)];
    var(i)=randi([1 Var],1);
end

for i=1:time_instances
    for j=1:mobile_nodes
        positions(j,:,i+1)=positions(j,:,i)+Delta*var(j)*[randn randn]*F;
        if positions(j,1,i+1)>(xyDim) 
           positions(j,:,i+1) =[positions(j,1,i+1)-xyDim,positions(j,2,i+1)];
        elseif positions(j,1,i+1)<0
           positions(j,:,i+1) =[positions(j,1,i+1)+xyDim,positions(j,2,i+1)];
        end
        if positions(j,2,i+1)>(xyDim)
            positions(j,:,i+1) =[positions(j,1,i+1),positions(j,2,i+1)-xyDim];
        elseif positions(j,2,i+1)<0
            positions(j,:,i+1) =[positions(j,1,i+1),positions(j,2,i+1)+xyDim];
        end
    end
end


Nodes=zeros(length(Anchors)+mobile_nodes,2,time_instances+1);
for i=1:time_instances+1
    Nodes(:,:,i)=[positions(:,:,i); Anchors];
end

Distance=zeros(mobile_nodes,length(Anchors)+mobile_nodes,time_instances+1);
Probabilities=zeros(mobile_nodes,length(Anchors)+mobile_nodes,time_instances+1);

for L=1:time_instances+1
    k=1;
    for i=1:mobile_nodes
        k=k+1;
        for j=k:length(Anchors)+mobile_nodes
            Distance(i,j,L)=abs(norm(positions(i,:,L)-Nodes(j,:,L))+randn);
            Probabilities(i,j,L)=exp(-(Distance(i,j,L)^2)./(2*(R^2)));
            
        end
    end

end

end
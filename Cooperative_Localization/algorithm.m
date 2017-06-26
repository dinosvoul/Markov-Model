clc; clear all;
%%
R=10;
space=1;
Var=0.5;
xyDim=20;
time_instances=5;
mobile_nodes=2;


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

%% 
%Mobility model.
for i=1:mobile_nodes
    positions(i,:,1)=[randi([0 xyDim],1),randi([0 xyDim],1)];
    var(i)=1;%randi([1 Var],1);
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
%Stack the Nodes together;
for i=1:time_instances+1
    Nodes(:,:,i)=[positions(:,:,i); Anchors];
end
%% 
%Range model.
Distance=zeros(length(Anchors)+mobile_nodes,length(Anchors)+mobile_nodes,time_instances+1);
Probabilities=zeros(mobile_nodes,length(Anchors)+mobile_nodes,time_instances+1);
noise=0.1;
for L=1:time_instances+1
    k=1;
    for i=1:mobile_nodes
        k=k+1;
        for j=k:length(Anchors)+mobile_nodes
            Distance(i,j,L)=abs(norm(positions(i,:,L)-Nodes(j,:,L))+noise*randn);
            Probabilities(i,j,L)=exp(-(Distance(i,j,L)^2)./(2*(R^2)));
            
        end
    end

end
%% 
%Connectivity model.
min_dist=min(min(min(nonzeros(Distance))));
ratio=savic(R,min_dist)./(1-savic(R,min_dist));

if lt(1,ratio)
    radians=atan(ratio);
    alpha_max1=tan(pi./2-radians);
    alpha=alpha_max1;%*rand;
else
    alpha_max1=1;
    alpha=alpha_max1*rand;
end

Connectivity_status=zeros(size(Distance,2),size(Distance,2),time_instances+1);
Connectivity_status(1:mobile_nodes,:,1)=Probabilities(1:mobile_nodes,:,1)>rand*ones(mobile_nodes, size(Probabilities,2));
for L=2:time_instances+1     
    k=1;
    for i=1:mobile_nodes
        k=k+1;
        for j=k:size(Distance,2)       
            prop=Probabilities(i,j,L);   
            K=[1-alpha,alpha];
            beta=(prop*alpha)./(1-prop);          
            L1=[beta,1-beta];
            Transition=[K;L1];
            if Connectivity_status(i,j,L-1)==1
                Connectivity_status(i,j,L)=Transition(1,1)>rand;    
            else
                Connectivity_status(i,j,L)=Transition(2,1)>rand;
            end 
        end
       
    end    
%        figure('units','normalized','outerposition',[0 0 1 1])
%     scatter(Anchors(:,1),Anchors(:,2),'filled');hold on;
% %     labels_1 = num2str((1:length(Anchors))','%d');
%  %   text(Anchors(:,1), Anchors(:,1), labels_1, 'horizontal','left', 'vertical','bottom');hold on;
%     labels = num2str((1:mobile_nodes)','%d');
%     text(positions(:,1,L), positions(:,2,L), labels, 'horizontal','left', 'vertical','bottom');hold on;
%     scatter(positions(:,1,L),positions(:,2,L),'g','filled')
%     gplot(Connectivity_status(:,:,L),Nodes(:,:,L),'-k')
%     ylim([ min(Nodes(:,2,L))-1 max(Nodes(:,2,L))+1.1])
%     xlim([ min(Nodes(:,1,L))-1 max(Nodes(:,1,L))+1.1])
%     axis off
%      pause(1)
%      close
end


%%
% Localization
es_pos=zeros(mobile_nodes,2);
Connection=zeros(size(Distance,2),size(Distance,2),time_instances+1);
Measurements=zeros(size(Distance,2),size(Distance,2),time_instances+1);
for L=1:time_instances+1
 Connection(:,:,L)= Connectivity_status(:,:,L)+Connectivity_status(:,:,L)';
 Measurements(:,:,L)=Distance(:,:,L)+Distance(:,:,L)';
end
fig=figure;
close
for L=1
    for i=1:mobile_nodes
        Observe= nonzeros(Connection(i,:,L).*Measurements(i,:,L));
        P=find(Connection(i,:,L));
        K=Nodes(P,:,L);
        Z=[Nodes(P,:,1)';zeros(1,length(P))]; 
        minX1=min(K(:,1));
        minX2=min(K(:,2));
        maxX1=max(K(:,1));
        maxX2=max(K(:,2));
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        lb = [minX1,minX2];
        ub = [maxX1,maxX2];
        fun=@(x) sum((Observe-sqrt((K(:,1)-x(1)).^2+(K(:,2)-x(2)).^2)).^2);
        %[sol,fval] = fminunc(fun,N1(2:3));
        sol=positions(i,:,1)+var(i)*[randn randn];
        sol1 = fmincon(fun,sol,A,b,Aeq,beq,lb,ub);
        es_pos(i,:,L)=sol1;
    end
%      figure('units','normalized','outerposition',[0 0 1 1])
%      scatter(Anchors(:,1),Anchors(:,2),'filled');hold on;
%     labels = num2str((1:mobile_nodes)','%d');
%     text(positions(:,1,L), positions(:,2,L), labels, 'horizontal','left', 'vertical','bottom');hold on;
%     scatter(positions(:,1,L),positions(:,2,L),'g','filled')
%     gplot(Connectivity_status(:,:,L),Nodes(:,:,L),'-k')
%     ylim([ min(Nodes(:,2,L))-1 max(Nodes(:,2,L))+1.1])
%     xlim([ min(Nodes(:,1,L))-1 max(Nodes(:,1,L))+1.1])
%     labels = num2str((1:mobile_nodes)','%d');
%     text(es_pos(:,1,L), es_pos(:,2,L), labels, 'horizontal','left', 'vertical','bottom');hold on;
%     scatter(es_pos(:,1,L),es_pos(:,2,L),'x','LineWidth',2.5)%      scatter(N1(2),N1(3),'k','filled')
%     pause(Delta)
%      close
end
%%
 
 
%  writerObj = VideoWriter('out.avi'); % Name it.
% writerObj.FrameRate = 1; % How many frames per second.
% open(writerObj); 

 % How many frames per second.
 for L=2:time_instances+1
    for i=1:mobile_nodes
        Observe= nonzeros(Connection(i,:,L).*Measurements(i,:,L));
        P=Connection(i,:,L)==1;
        K=Nodes(P,:,L);
        minX1=min(K(:,1));
        minX2=min(K(:,2));
        maxX1=max(K(:,1));
        maxX2=max(K(:,2));
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        lb = [minX1,minX2];
        ub = [maxX1,maxX2];
        sol=es_pos(i,:,L-1)+Delta*[var(i) var(i)]*F;
        fun=@(x) sum((Observe-sqrt((K(:,1)-x(1)).^2+(K(:,2)-x(2)).^2)).^2);
        sol1 = fmincon(fun,sol,A,b,Aeq,beq,lb,ub);
        es_pos(i,:,L)=sol1;
    end
%     clf(fig);
%      figure('units','normalized','outerposition',[0 0 1 1])
   figure(L-1)
%     hold on;
     
        
      scatter(Anchors(:,1),Anchors(:,2),'filled');
     labels = num2str((1:mobile_nodes)','%d');
     text(positions(:,1,L), positions(:,2,L), labels, 'horizontal','left', 'vertical','bottom');hold on;
    scatter(positions(:,1,L),positions(:,2,L),'g','filled')
%     gplot(Connectivity_status(:,:,L),Nodes(:,:,L),'-k')
     ylim([ min(Nodes(:,2,L))-1 max(Nodes(:,2,L))+1.1])
     xlim([ min(Nodes(:,1,L))-1 max(Nodes(:,1,L))+1.1])
     labels = num2str((1:mobile_nodes)','%d');
%      text(es_pos(:,1,L), es_pos(:,2,L), labels, 'horizontal','left', 'vertical','bottom');hold on;
    scatter(es_pos(:,1,L),es_pos(:,2,L),'r','LineWidth',2.5)%      scatter(N1(2),N1(3),'k','filled')
%      hold off
        axis off
     pause(1)
%     frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%         writeVideo(writerObj, frame);
    
     close    
%     F=getframe(gcf);
 end
%hold off
% close(writerObj); % Saves the movie.

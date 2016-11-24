function [Connectivity_status]=connection(Distance,probabilities,time_instances,alpha,mobile_nodes)
Connectivity_status=zeros(mobile_nodes,size(Distance,2),time_instances+1);
Connectivity_status(1:mobile_nodes,:,1)=probabilities(1:mobile_nodes,:,1)>rand*ones(mobile_nodes, size(probabilities,2));
for L=2:time_instances+1
    k=1;
    for i=1:mobile_nodes
        k=k+1;
        for j=k:size(Distance,2)       
            prop=probabilities(i,j,L);   
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
end
end
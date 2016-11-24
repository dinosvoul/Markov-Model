function [Result]=factorize(R,alpha_space,Markov,Distance)
    min_distance=min(min(Distance));
    frac=exp(-min_distance^2/(2*R^2))./(1-exp(-min_distance^2/(2*R^2)));
    
    if lt(1,frac)
        radians=atan(frac);
        alpha_max=tan(pi./2-radians);        
    else
        alpha_max=1;        
    end
    synolo=zeros(max(size(alpha_space)),max(size(Markov)));
    if Markov(1)==1
        synolo(:,1)=exp(-Distance(1)^2/(2*R^2))*ones(max(size(alpha_space)),1); 
        
    else
        synolo(:,1)=(1-exp(-Distance(1)^2/(2*R^2)))*ones(max(size(alpha_space)),1);
       
    end
    
    for i=2:max(size(Markov))
         rat=exp(-Distance(i)^2/(2*R^2))./(1-exp(-Distance(i)^2/(2*R^2)));
         if Markov(i)==1 && Markov(i-1)==0         
            synolo(:,i)=alpha_space.*(rat); 
            
        else if Markov(i)==1 && Markov(i-1)==1
            synolo(:,i)=ones(max(size(alpha_space)),1)-(alpha_space); 
            
        else if Markov(i)==0 && Markov(i-1)==1
                synolo(:,i)=alpha_space;  
               
        else if Markov(i)==0 && Markov(i-1)==0
                synolo(:,i)=ones(max(size(alpha_space)),1)-(alpha_space.*rat);
                
             end
             end
             end
        end
    end
   in=find(alpha_space>alpha_max);
   
   if numel(in);not(0)
    synolo(in(1):max(size(alpha_space)),:)=ones(max(size(alpha_space))-in(1)+1,max(size(Markov)));
   end
   Result=synolo;    
end
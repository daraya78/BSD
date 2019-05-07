 function miplot(type,data,sreal,sest)
         
            T=size(data,1);
            ndim=size(data,2);
            state_un=unique(sreal);
            state_dim=size(state_un,2);
                  
           %Color definition
            colors=[[0 0.5 1];[0 1 1];[1 1 0];[1 0 0];[0 0.5 1];[0.5 1 0.5];[1 0.5 0]];
            colorcur=['r' 'b' 'g'];
            colormap(colors); 
           %----------------------
            % plot state sequence
            ma=max(max(data));
            mi=min(min(data));
            subplot(9,1,[1:4]);
            h1=image(sreal);
            h1.Parent.YLim=[mi ma];
            h1.YData=[mi ma];
            hold on
            plot(data(:,1),colorcur(1))
            grid
            plot(data(:,2),colorcur(2))
            grid
            title('Secuencia de Estados Original')
            xlabel('tiempo ms')
            ylabel('Señal')
            
            subplot(9,1,[6:9]);
            h2=image(sest);
            h2.Parent.YLim=[mi ma];
            h2.YData=[mi ma];
            hold on
            plot(data(:,1),colorcur(1))
            grid
            plot(data(:,2),colorcur(2))
            grid
            title('Secuencia de Estados Estimada')
            xlabel('tiempo ms')
            ylabel('Señal')
            
            legend('Fuente 1','Fuente 2')
            

   
 end
           
     
    
clear all   
% X_0 = [7,9.4,12.5,14,15.9,19.3,24.1,25.8,28.7,39.6,42.2,58.3,77.5,89.6,98,106.4];
   % X_0 = [0,1931,1724,1517,1345,1207,1069,952,848,745,669];
   %X_0 = [1320,1881,2438,2664,2754,2934,2949,2560,2452,2287];
    %X_0 = [1881,2438,2664,2754,2934,2949,2560,2452,2287];
    X_0=[106859,110980,127159,137172,148900,165949,212474,237704,252068,273600,301229,319067,333405,350700,355823,0,0,0];

%sol = [1,9,0.999902421,1.010715057];
sol = [1,9,0.72484,1.08363];

        temp=0;
        RGM_x_0=X_0;
        RGM_x_1=X_0;
        X_P=X_0;
        for i=sol(1):sol(2)
             temp=temp+RGM_x_0(i);
             RGM_x_1(i)=temp;
        end
        
        k=1;
        for i=sol(1):sol(2)-1
            %RGM_B(k,1)=-1/2.*(RGM_x_1(i)+RGM_x_1(i+1));
            RGM_B(k,1)=-1*((sol(3))*RGM_x_1(i)+(1-sol(3))*RGM_x_1(i+1));
          k=k+1;
        end
        
        RGM_B(:,2)=sol(4);
        RGM_y=(RGM_x_0(sol(1)+1:sol(2)))';
        RGM_v=inv(RGM_B'*RGM_B)*RGM_B'*RGM_y;
        for i=sol(1)+1:sol(2)+2
            X_P(i)=(X_0(sol(1))-(RGM_v(2)/RGM_v(1)))*exp(-(RGM_v(1)*(i-1)))*(1-exp(RGM_v(1)));
        end
        
        X_0 = [2080,2261];
        X_P1 = X_P(10:11);
        fprintf('%7.0f  ',X_P);
        fprintf('\n');
        fprintf('E1 = %10.5f  \n',mape(X_0, X_P1));
        fprintf('%7.0f  ',X_P1);

 function Ans = mape( Y, Ypredict)
smape = 0;
        for i = 1 :length(Y)
        if (Y(i)~=0)
            smape = smape + (abs((Ypredict(i) - Y(i))) / Y(i));
        end
        end
Ans = smape * 100/length(Y);
 end
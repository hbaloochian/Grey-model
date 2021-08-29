function [BestGen,BestFitness,gx]=HarmonySearch

global NVAR NG NH MaxItr HMS HMCR PARmin PARmax bwmin bwmax;
global HM NCHV fitness PVB BW gx;
global BestIndex WorstIndex BestFit WorstFit BestGen currentIteration;
global X_0 SNUMBER IT PER; 


NVAR=4;         %number of variables
NG=6;           %number of ineguality constraints
NH=0;           %number of eguality constraints
MaxItr=5000;    % maximum number of iterations
HMS=50;          % harmony memory size
HMCR=0.7;       % harmony consideration rate  0< HMCR <1
PARmin=0.4;      % minumum pitch adjusting rate
PARmax=0.9;      % maximum pitch adjusting rate
bwmin=0.0001;    % minumum bandwidth
bwmax=1.0;      % maxiumum bandwidth
SNUMBER=9;

PVB=[1 SNUMBER;1 SNUMBER;0.0 1;0.0 2];   % range of variables
IT=0;

PER=2;  %number of Prediction

%********************************************************
%                  Add your series here  
X_0 = [1881,2438,2664,2754,2934,2949,2560,2452,2287];
%********************************************************


% /**** Initiate Matrix ****/
HM=zeros(HMS,NVAR);
NCHV=zeros(1,NVAR);
BestGen=zeros(1,NVAR);
fitness=zeros(1,HMS);
BW=zeros(1,NVAR);
gx=zeros(1,NG);


% warning off MATLAB:m_warning_end_without_block

MainHarmony;

% /**********************************************/
    function sum =FGM_1_1(sol)
    %sol = [1,6,0.5,1];
        
        X_0(1) = 0;
        temp=0;
        FGM_x_0=X_0;
        FGM_x_1=X_0;
        X_P=X_0;
        for i=sol(1):sol(2)
             temp=temp+FGM_x_0(i);
             FGM_x_1(i)=temp;
        end
        
        k=1;
        for i=sol(1):sol(2)-1
            FGM_B(k,1)=-1*((sol(3))*FGM_x_1(i)+(1-sol(3))*FGM_x_1(i+1));
          k=k+1;
        end
        
        FGM_B(:,2)=sol(4);
        FGM_y=(FGM_x_0(sol(1)+1:sol(2)))';
        FGM_v=inv(FGM_B'*FGM_B)*FGM_B'*FGM_y;
        for i=sol(1)+1:sol(2)
            X_P(i)=(X_0(sol(1))-(FGM_v(2)/FGM_v(1)))*exp(-(FGM_v(1)*(i-1)))*(1-exp(FGM_v(1)));
        end
        fprintf('E1 = %10.5f  \n',mape(X_0, X_P));
        fprintf('%7.5f  ',X_P);
        sum = mape(X_0, X_P);
    end

% /**********************************************/
function Ans = mape( Y, Ypredict)
smape = 0;
        for i = 1 :length(Y)
            if (Y(i) ~= 0)
                smape = smape + (abs((Ypredict(i) - Y(i))) / Y(i));
            end
        end
Ans = smape * 100/length(Y);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function sum =Fitness(sol)
        
        X_0(1) = 0;
        temp=0;
        FGM_x_0=X_0;
        FGM_x_1=X_0;
        X_P=X_0;
        for i=sol(1):sol(2)
             temp=temp+FGM_x_0(i);
             FGM_x_1(i)=temp;
        end
        
        k=1;
        for i=sol(1):sol(2)-1
            FGM_B(k,1)=-1*((sol(3))*FGM_x_1(i)+(1-sol(3))*FGM_x_1(i+1));
          k=k+1;
        end
        
        FGM_B(:,2)=sol(4);
        FGM_y=(FGM_x_0(sol(1)+1:sol(2)))';
        FGM_v=inv(FGM_B'*FGM_B)*FGM_B'*FGM_y;
        for i=sol(1)+1:sol(2)
            X_P(i)=(X_0(sol(1))-(FGM_v(2)/FGM_v(1)))*exp(-(FGM_v(1)*(i-1)))*(1-exp(FGM_v(1)));
        end

        sum = mape(X_0, X_P);

     end

% /*********************************************/
    function sum=eg(sol)
        
        % constraints g(x) > 0
        gx(1)=sol(1)-0.0193*sol(3);     %  x1 - 0.0193 x3 > 0
        gx(2)=sol(2)-0.00954*sol(3);
        gx(3)=3.14*sol(3)^2*sol(4)+(4/3)*3.14*sol(3)^3 - 1296000;
        gx(4)=-sol(4)+240;
        gx(5)=sol(1) - 1.1;
        gx(6)=sol(2) - 0.6;
        
        % we use static penalty function to handle constraints
        sum = 0;
        for i=1:NG
            if(gx(i)<0)
                sum = sum - 1000 * gx(i);
            end
        end
    end

% /*********************************************/

    function initialize
        % randomly initialize the HM
        for i=1:HMS
            HM(i,1)=1;
            HM(i,2)=SNUMBER;
            HM(i,3) = randval( PVB(3,1), PVB(3,2));  %rand(1);
            HM(i,4) = randval( PVB(4,1), PVB(4,2));  %rand(1);
            fitness(i) = Fitness(HM(i,:));
        end
    end

%/*******************************************/

    function MainHarmony
        % global NVAR NG NH MaxItr HMS HMCR PARmin PARmax bwmin bwmax;
        % global HM NCHV fitness PVB BW gx currentIteration;
        
        initialize;
        currentIteration  = 0;
        
        while(StopCondition(currentIteration))
            
            PAR=(PARmax-PARmin)/(MaxItr)*currentIteration+PARmin;
            coef=log(bwmin/bwmax)/MaxItr;
            for pp =1:NVAR
                BW(pp)=bwmax*exp(coef*currentIteration);
            end
            % improvise a new harmony vector
            for i =3:NVAR
                ran = rand(1);
                if( ran < HMCR ) % memory consideration
                    index = randint(1,HMS);
                    NCHV(i) = HM(index,i);
                    pvbRan = rand(1);
                    if( pvbRan < PAR) % pitch adjusting
                        pvbRan1 = rand(1);
                        result = NCHV(i);
                        if( pvbRan1 < 0.5)
                            result = result+  rand(1) * BW(i);
                            if( result < PVB(i,2))
                                NCHV(i) = result;
                            end
                        else
                            result = result- rand(1) * BW(i);
                            if( result > PVB(i,1))
                                NCHV(i) = result;
                            end
                        end
                    end
                else
                    NCHV(i) = randval( PVB(i,1), PVB(i,2) ); % random selection
                end
            end
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %  NCHV = [1,SNUMBER,0.5,1];     
            NCHV(1) = 1;     
            NCHV(2) = SNUMBER;     
            newFitness = Fitness(NCHV);
            UpdateHM( newFitness );
            currentIteration=currentIteration+1;
        end
        BestFitness = min(fitness);
        XX=HM;
        for i=1:size(fitness,2)
            if (fitness(i)==BestFitness)
                BestIndex = i;
            end
            XX(i,5)=fitness(i);
            XX(i,6)=i;
        end
        fprintf('-------------------Ansser Is--------------\n');
        fprintf('I   J    Alpha  Fitness Index\n');
        fprintf('%d   %d   %7.5f    %7.5f     %7.5f   %d\n',XX');
        fprintf('\n');
        fprintf('------------------------------------------\n');
        fprintf('BestFitness = %7.5f',BestFitness);
        fprintf('\n'); 
        fprintf('BestIndex = %d',BestIndex);
        fprintf('\n');
        fprintf('SERI =');
        fprintf('  %7.5f',X_0(HM(BestIndex,1):HM(BestIndex,2)));
        fprintf('\n'); 
        fprintf('Index = %d   %d     Alpha=%7.5f   Beta=%7.5f',HM(BestIndex,1),HM(BestIndex,2),HM(BestIndex,3),HM(BestIndex,4));
        fprintf('\n'); 
        v=[1,SNUMBER,0.5,1];
        FGM_1_1(v);
        fprintf(' Orginal \n'); 
        v=[1,SNUMBER,HM(BestIndex,3),HM(BestIndex,4)];
        FGM_1_1(v);
        fprintf(' Peridict \n'); 
        %v=[1,SNUMBER,0.44146,0.95287];
        %GM_1_1(v);
    end
% /*****************************************/

    function UpdateHM( NewFit )
        % global NVAR MaxItr HMS ;
        % global HM NCHV BestGen fitness ;
        % global BestIndex WorstIndex BestFit WorstFit currentIteration;
        
        if(currentIteration==0)
            BestFit=fitness(1);
            for i = 1:HMS
                if( fitness(i) < BestFit )
                    BestFit = fitness(i);
                    BestIndex =i;
                    IT = currentIteration;
                end
            end
            
            WorstFit=fitness(1);
            for i = 1:HMS
                if( fitness(i) > WorstFit )
                    WorstFit = fitness(i);
                    WorstIndex =i;
                end
            end
        end
        if (NewFit< WorstFit)
            
            if( NewFit < BestFit )
                HM(WorstIndex,:)=NCHV;
                BestGen=NCHV;
                fitness(WorstIndex)=NewFit;
                BestIndex=WorstIndex;
                %IT = currentIteration;
            else
                HM(WorstIndex,:)=NCHV;
                fitness(WorstIndex)=NewFit;
                %IT = currentIteration;
            end
            
            
            WorstFit=fitness(1);
            WorstIndex =1;
            for i = 1:HMS
                if( fitness(i) > WorstFit )
                    WorstFit = fitness(i);
                    WorstIndex =i;
                end
            end
            
        end
    end % main if
end %function

% /*****************************************/
function val1=randval(Maxv,Minv)
    val1=rand(1)*(Maxv-Minv)+Minv;
end

function val2=randint(Maxv,Minv)
    val2=round(rand(1)*(Maxv-Minv)+Minv);
end
% /*******************************************/

function val=StopCondition(Itr)
    global MaxItr;
    val=1;
    if(Itr>MaxItr)
        val=0;
    end
end

% /*******************************************/





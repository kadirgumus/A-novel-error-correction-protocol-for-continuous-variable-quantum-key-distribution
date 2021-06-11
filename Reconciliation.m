%% Defining channel parameters
EsN0dB = -15.23; %SNR
EsN0 = 10^(EsN0dB/10);
T = (EsN0*(1+0.1)/0.5)/(0.5-0.01*EsN0);
C = 0.5*log2(1+EsN0);
%% Defining parameters for the simulation
Simulation_Iterations = 10;

%% Defining code parameters for the base code
N = 10^6; %Block length
p = 0; %Amount of punctured bits, increases the rate of the code
s = 0; %Amount of shortened bits, decreases the rate of the code
Rate = 0.02; %The rate of the base code, currently only 0.02 and 0.1 are available

%% Defining parameters for the new reconciliation
Decoding_Attempts_Max = 2; %The maximum amount of decoding  
s_protocol = 0.06*N*Rate; %The amount of bits revealed during each decoding attempt, an array of length decoding attempts max -1
p_protocol = 0; %The amount of bits shortened during eacht decoding attempt, should always be 0
Beta = [Rate (Rate - s_protocol/N)]/C;
%% Defining the parameters for the decoding
MaxIter = [500; 500]; %The maximum amount of decoding iterations for the decoding attempts, should be a matrix of size Decoding_attempts_max x n, where n is an arbitrary number and corrsponds to the amount of decoding iteration intervals
parameters.OutputValue= 'Information Part'; %'Information Part' for just the information bits, 'Whole Codeword' for the entire codeword
parameters.DecisionMethod = 'Hard Decision'; %'Hard Decision' for bits as output, 'Soft Decision' for LLRs as output
parameters.maxUnchangedIterations = 5; %0 If you do not want early termination based on stagnating LLRs, any positive number X if you want the decoding to step after X decoding iterations without change in output bits
parameters.Algorithm = 'Flooding'; %'Flooding' for flooding decoder and 'Layered' for layered decoding
parameters.Algorithm2 = 'Sum Product'; %'Sum Product' for the sum product algorithm, 'Min Sum' for the min sum algoritm
parameters.Threads = 1; %The amount of threads you want to use for the multithread calculation, only relevant for the Flooding algorithm
parameters.normfactor= 1 ; %The normalization factor for normalization min-sum decoding. 1 if you want normal min sum <1 if you want offset decoding
parameters.N = N;
parameters.Rate = Rate;
%% Allocating memory for variables
errors = zeros(size(MaxIter,2),Simulation_Iterations,Decoding_Attempts_Max);
Decoding_Iterations = zeros(size(MaxIter,2),Simulation_Iterations,Decoding_Attempts_Max);
FER = zeros(size(MaxIter,2),Decoding_Attempts_Max);
SKR = zeros(size(MaxIter,2),Decoding_Attempts_Max);
%% Running the simulations
for i = 1:Simulation_Iterations
    %Generating Alice's and Bob's qunatum states
    X = sqrt(2)*0.5*(2*round(rand([N,1]))-1); %Alice's generated quantum states
    EsN0_1D = EsN0/mean(X.^2);
    EsN0dB_1D = 10*log(EsN0_1D)/log(10);
    Z = normrnd(0,1/sqrt(EsN0_1D),[N 1]); %Quantum channel noise
    Y = X+Z; %Bob's received quantum states
    
    %Generating a random bit sequence
    S = logical(randi([0 1], Rate*N, 1)); 
    
    %Generating a random parity check matrix
    if mod(i,100) == 1 
        if Rate == 0.02
            [Hpos,H] = LDPC_Generator(N,[9/400 7/400],[2 3; 57 57],[17/1600 3/320 3/5 9/25],[3 7 0 0 ;0 0 2 3],'T',0.02);
        elseif Rate == 0.1
            [Hpos,H] = LDPC_Generator(N,[0.1 0.025],[3 3; 20 25],[0.025 0.875],[15 0;0 3],'T',0.1);
        end
       hEnc = comm.LDPCEncoder(H);
    end
    %Encoding the codeword 
    C = step(hEnc,S); 
    
    %Generating which bits are punctured and shortened
    p_s_idx = randperm(N*Rate,p+s+sum(s_protocol)+sum(p_protocol)); 
    p_idx = p_s_idx(1:p);
    s_idx = p_s_idx(p+1:p+s);
    Shortened_bits = C(s_idx); 
    for d = 1:Decoding_Attempts_Max - 1
        p_idx_protocol{d} = p_s_idx(p+s+sum(p_protocol(1:d-1))+sum(s_protocol(1:d-1))+1:p+s+sum(p_protocol(1:d))+sum(s_protocol(1:d-1)));
        s_idx_protocol{d} = p_s_idx(p+s+sum(p_protocol(1:d))+sum(s_protocol(1:d-1))+1:p+s+sum(p_protocol(1:d))+sum(s_protocol(1:d)));
        Shortened_bits_protocol{d} = C(s_idx_protocol{d});
    end
    %Determining the LLRs
    l = Multi_Dimensional_Reconciliation(1,C,X,Y,EsN0dB_1D,1); 
    l = LLRps(l,p_idx,s_idx,Shortened_bits);
    Lq_in = [];
    Lq_out2 = [];
    for j = 1:size(MaxIter,2)
        if strcmp(parameters.Algorithm,'Layered') == 1 && j > 1
            l_temp = l_temp2;
        else
            l_temp = l;
        end
        hError = comm.ErrorRate;
        for c = 1:Decoding_Attempts_Max
            if c == 1
                if j == 1
                    parameters.MaximumIterationCount = MaxIter(c,j);
                else
                    parameters.MaximumIterationCount = MaxIter(c,j)-MaxIter(c,j-1);
                end
            else
                parameters.MaximumIterationCount = MaxIter(c,j);
            end
            [S_hat, iterations, Lq_out, parityCheck, LLR_out] = LDPC_Decoder(l_temp, Hpos, parameters, Lq_in);          
            errorStats     = step(hError, S, logical(S_hat'));
            fprintf('Error rate       = %1.2f\nNumber of errors = %d\n',errorStats(1), errorStats(2))
            errors(j,i,c) = errorStats(2);
            Decoding_Iterations(j,i,c) = iterations;
            if  strcmp(parameters.Algorithm,'Layered') == 1
                l_temp = LLR_out;
                l_temp2 = l_temp;
            end
            if c == 1
                Lq_temp = Lq_out;
            end
            if nnz(parityCheck)~= 0 && errorStats(2) == 0 
                errors(j,i,c) = -1;
            elseif nnz(parityCheck)== 0 && errorStats(2) == 0 
                break;
            end
            if Decoding_Attempts_Max > 1 && c < Decoding_Attempts_Max
                hError = comm.ErrorRate;
                l_temp = LLRps(l_temp,p_idx_protocol{c},s_idx_protocol{c},Shortened_bits_protocol{c});
            end
            Lq_in = Lq_out;
         end
        Lq_in = Lq_temp;
    end
end

for i = 1:size(MaxIter,2)
    for j = 1:Decoding_Attempts_Max
        FER(i,j) = nnz(errors(i,:,j))/Simulation_Iterations;
        if j == 1
            [~,SKR(i,j),~] = calcSKR( 0.5, 0.5, T, 0.01, 0.1, Beta(j), FER(i,j));
        else
            [~,temp,~] = calcSKR( 0.5, 0.5, T, 0.01, 0.1, Beta(j), FER(i,j));
            SKR(i,j) = SKR(i,j-1) + FER(i,j-1)*temp;
        end
    end
end
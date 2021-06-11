function [H_pos, H] = LDPC_Generator(N,VN_dist,VN_amount, CN_dist, CN_amount,Type,Rate,q)
%N: The length of the codeword
%VN_dist: An array conataining the distributions probabilities of the different VN
%In the case of Type 'R' or 'QS' do not include the nodes in the diagional
%part with only one connection.
%VN_amount: A matrix containing the amount of edges a VN contains per edge type
%CN_dist: An array conataining the distributions probabilities of the different CN
%CN_amount:  A matrix containing the amount of edges a CN contains per edge type
%In the case of type 'R' or 'QS' do not include the edge connections to
%VN's with degree 1. 
%Type: There are 3: 'T' = Triangular, with a structure described by Urbanke, 'QS' = Quasi-Cyclic (Milicevic 2017), 
%'R' = random, so without a lower triangular structure, 
%not adviced as this makes encoding much more difficult. 
%Rate: The Rate of your code, not required for 'TR'
%q: only relevant for QS codes
if Type == 'R'
    H_pos = zeros(round(sum(sum(VN_amount).*VN_dist)*N),2);
    Omega = zeros(N,size(VN_amount,1));
    idx = 0;
    for i = 1:length(VN_dist)
        for j = 1:size(VN_amount,1)
            Omega(N*idx+1:N*(idx+VN_dist(i)),j) = VN_amount(j,i)*ones(VN_dist(i)*N,1);
        end
        idx = idx + VN_dist(i);
    end
    for i = 1:size(CN_amount,2)
        temp = [];
        idx = 1;
        for j = 1:length(CN_dist)
            temp = [temp repmat(idx:idx+CN_dist(j)*N-1,[1 CN_amount(j,i)])];
            idx = idx + N*CN_dist(j);
        end
        possible_idx{i} = temp(randperm(length(temp)));
    end
    sums = ones([size(VN_amount,1) 1]);
    idx_H = 1;
    for j = 1:N
        check = 0;
        check2 = 0;
        while check ==0
            if check2 == 500
                error('The algorithm could not find an LDPC code this time, try again.');
            end
            H_idx = [];
            for i = 1:size(VN_amount,2)
                idx = possible_idx{i}(sums(i):sums(i)+Omega(j,i)-1);
                H_idx = [H_idx idx];
            end
            if length(unique(H_idx)) == sum(Omega(j,:))
                check = 1;
            else
                for i = 1:size(VN_amount,2)
                    temp = possible_idx{i}(sums(i):end);
                    possible_idx{i}(sums(i):end) = temp(randperm(length(temp)));
                end
            end
            check2 = check2+1;
        end
        H_pos(idx_H:idx_H+length(H_idx)-1,:) = [H_idx;j*ones(length(H_idx),1)']';
        sums = sums + Omega(j,:)';
        idx_H = idx_H + sum(Omega(j,:));  
    end
    H = sparse(H_pos(:,1),H_pos(:,2),ones(length(H_pos),1));
elseif Type == 'T'
    check3 = 0;
    H_pos = zeros(round(sum(sum(VN_amount).*VN_dist)*N) + round((1-sum(VN_dist))*N) - VN_amount(1,2)*(VN_amount(1,2)-1),2);
    H_pos(1:round((1-Rate)*N),:) = [1:round((1-Rate)*N); round(Rate*N)+1:N]';
    Omega = zeros(round(sum(VN_dist)*N),size(VN_amount,1));
    idx = 0;
    for i = 1:length(VN_dist)
        for j = 1:size(VN_amount,1)
            Omega(round(N*idx)+1:round(N*(idx+VN_dist(i))),j) = VN_amount(j,i)*ones(round(VN_dist(i)*N),1);
        end
        idx = idx + VN_dist(i);
    end
    Omega(round(sum(VN_dist)*N)-VN_amount(1,2)*(VN_amount(1,2)-1)+1:end,1) = 0;
    Omega(round(Rate*N)+1:round(sum(VN_dist)*N)-VN_amount(1,2)*(VN_amount(1,2)-1),1) = Omega(round(Rate*N)+1:round(sum(VN_dist)*N)-VN_amount(1,2)*(VN_amount(1,2)-1),1)-1;
    Use_Check = [];
    for i = 1:nnz(CN_amount(1,:))
        Use_Check = [Use_Check (CN_amount(1,i)-1)*ones(1,round(CN_dist(i)*N))];
    end

    idx_H = round((1-Rate)*N)+1;
    for j = fliplr(round(Rate*N):round(sum(VN_dist)*N)-VN_amount(1,2)*(VN_amount(1,2)-1))
        range = j-round(N*Rate)+1:round(((sum(VN_dist)-Rate)*N));
        temp_check = Use_Check( j-round(N*Rate)+1:round(((sum(VN_dist)-Rate)*N)));
        idx = randsample(range(temp_check ~= 0),Omega(j,1));
        Use_Check(idx) = Use_Check(idx)-1;
        H_pos(idx_H:idx_H+Omega(j,1)-1,:) = [idx' j*ones(Omega(j,1),1)];
        idx_H = idx_H + Omega(j,1);
        Omega(j,1) = 0;
    end
    temp = [];
    idx = 1;
    for i = 1:length(Use_Check)
        temp = [temp i*ones(1,Use_Check(i))];
    end

    possible_idx{1} = temp(randperm(length(temp)));
   
    for i = 2:size(CN_amount,1)
        temp = [];
        idx = 1;
        for j = 1:length(CN_dist)
            temp = [temp repmat(idx:idx+round(CN_dist(j)*N)-1,[1 CN_amount(i,j)])];
            idx = idx + round(N*CN_dist(j));
        end
        possible_idx{i} = temp(randperm(length(temp)));
    end
    sums = ones([size(VN_amount,2) 1]);
    for j = fliplr(1:length(Omega))
        check = 0;
        check2 = 0;
        while check ==0 && check3 == 0 
            if check2 == 500
                %error('The algorithm could not find an LDPC code this time, try again.');
                check3 = 1;
            end
            H_idx = [];
            for i = 1:size(VN_amount,2)
                idx = possible_idx{i}(sums(i):sums(i)+Omega(j,i)-1);
                H_idx = [H_idx idx];
            end
            if length(unique(H_idx)) == sum(Omega(j,:))
                 check = 1;
            else
                for i = 1:size(VN_amount,2)
                    temp = possible_idx{i}(sums(i):end);
                    possible_idx{i}(sums(i):end) = temp(randperm(length(temp)));
                end
            end
            check2 = check2+1;
        end
        H_pos(idx_H:idx_H+length(H_idx)-1,:) = [H_idx;j*ones(length(H_idx),1)']';
        sums = sums + Omega(j,:)';
        idx_H = idx_H + sum(Omega(j,:));  
    end
    H_pos = H_pos(H_pos(:,1) ~= 0,:);
    H_pos = H_pos(H_pos(:,1) <= floor((1-Rate)*N),:);
    H = sparse(H_pos(:,1),H_pos(:,2),ones(length(H_pos),1));
    if check3 == 1
        [H_pos, H] = LDPC_Generator(N,VN_dist,VN_amount, CN_dist, CN_amount,Type,Rate);
    end
elseif Type == 'QS'
    if mod(N/q,1)~=0
        error('N divided by q has to be an integer number');
    end
    [H_pos, ~] = LDPC_Generator(N/q,VN_dist,VN_amount,CN_dist,CN_amount,'T',Rate);
    I_num = [ones(round((1-Rate)*N/q),1);randi(q,[length(H_pos)-round((1-Rate)*N/q) 1])];
    H_pos_new = zeros([length(I_num)*q, 2]);
    for i = 1:length(I_num)
        H_pos_new(1+(i-1)*q:i*q,:) = repmat((H_pos(i,:)-1)*q,[q 1]) + [1:q;circshift(1:q,I_num(i)-1)]';
    end
    H = sparse(H_pos_new(:,1),H_pos_new(:,2),1);
    H_pos = H_pos_new;
end
return

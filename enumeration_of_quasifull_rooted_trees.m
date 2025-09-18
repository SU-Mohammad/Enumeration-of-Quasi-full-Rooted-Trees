%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Last modified on the September, 2025
%==========================================================================           
%          Completed on the October, 2023
% Copy right: S.U. Mohammad, Dept. of Mathematics, SUST 
%             Sylhet-3114, Bangladesh 
%==========================================================================
function ret = enumeration_of_quasifull_rooted_trees
%==========================================================================
% Returns the number of quasi-full t-ary rooted tress with number of nodes 
% (elements) num_of_nodes for all 2<=t<num_of_nodes. Note that here 
% num_of_nodes>=3 and ret is an array of num_of_nodes-2 elements.
num_of_nodes = input('Enter the number of elements (nodes): ');
ret = zeros(num_of_nodes-2,1);
disp(['Arity: ' 'Number of quasi-full rooted trees with ' ...
    num2str(num_of_nodes) ' elements']);
for arity = 2:num_of_nodes-1
    ret(arity-1) = count_quasifull_tary_rooted_trees(arity,num_of_nodes);
    disp([num2str(arity) ': ' num2str(ret(arity-1))]);
end
ret = sum(ret);
end
%==========================================================================
function ret = count_quasifull_tary_rooted_trees(arity,num_of_nodes)
%==========================================================================
% Returns the number of quasi-full rooted tress with arity t and number 
% of elements num_of_nodes. 
if arity+1==num_of_nodes
    ret = 1;
    return;
end
% Initializing our aassumption for the numbers of quasi-full t-ary ...
% rooted trees up to minimum number of elements t+1.
num_of_obj = ones(1,arity+1);
for n=arity+1:num_of_nodes-1
    ndid_lengths = floor(n/2);
    for m=1:arity-1
        ndid_lengths = nondec_inter_distant_lengths(ndid_lengths,n);
        % Due to problem with n=3 only.
        if ndid_lengths(1)==1&&m==n-1
            ndid_lengths(2) = 2;
        end
        t_ary_qtree = 0;
        for j=1:size(ndid_lengths,1) % value of p^
            lengths = [0 ndid_lengths(j,:) n];
            i = 2;
            prod = 1;
            while i<=m+2
                r_mjk = lengths(i)-lengths(i-1);
                t_mjk = 0; % repetition of same ordered terms.
                while i<=m+1 && r_mjk==lengths(i+1)-lengths(i)
                    t_mjk = t_mjk+1;
                    i = i+1;
                end
                prod = prod*nchoosek(num_of_obj(r_mjk)+t_mjk,1+t_mjk);
                i = i+1;
            end
            t_ary_qtree = t_ary_qtree+prod;
        end
        ret = t_ary_qtree;
    end
    num_of_obj = [num_of_obj ret];    
end
end
% =========================================================================
function ret = nondec_inter_distant_lengths(ndid_lengths,n)
%==========================================================================
if length(ndid_lengths)==1
    ret = (1:ndid_lengths)';
    return;
end
rownum = size(ndid_lengths,1);
colnum = size(ndid_lengths,2);
count = 1;
for i=1:rownum
    if colnum==1
        lengths = ndid_lengths(i);
        diff = lengths(1);
    else
        lengths = ndid_lengths(i,:);
        diff = lengths(colnum)-lengths(colnum-1);
    end
    flag = 1;
    while flag==1
        if lengths(colnum)+2*diff<=n
            lengths(colnum+1) = lengths(colnum)+diff;
            diff = diff+1;
            ret(count,:) = lengths;
            count = count+1;
        else
            flag = 0;
        end
    end
end
end

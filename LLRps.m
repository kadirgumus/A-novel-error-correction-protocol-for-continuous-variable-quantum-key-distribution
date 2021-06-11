function l = LLRps(l_org,p_idx,s_idx,Shortened_bits)
LLR_Shortened_bits = 10^10*(-1).^(Shortened_bits);
N = length(l_org);
l_idx = 1:N;
l_idx([p_idx s_idx]) = [];
l = 10^-7*ones([N 1]);
l(s_idx) = LLR_Shortened_bits;
l(l_idx) = l_org(l_idx);
end
function c = Octonion_mult(a,b)
%A function for multiplying octonions.
%Both a and b are matrices of size (N,8)
%c is the output octonion of size(N,8)
p = a(:,1:4);
q = a(:,5:8);
r = b(:,1:4);
s = b(:,5:8);
r_c = [r(:,1) -r(:,2) -r(:,3) -r(:,4)];
s_c = [s(:,1) -s(:,2) -s(:,3) -s(:,4)];
c1 = Quaternion_mult(p,r) - Quaternion_mult(s_c,q);
c2 = Quaternion_mult(s,p) + Quaternion_mult(q,r_c);
c = [c1 c2];
end
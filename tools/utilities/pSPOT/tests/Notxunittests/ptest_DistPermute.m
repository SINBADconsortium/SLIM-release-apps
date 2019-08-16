function test_DistPermute

x = parpool_size();
assert(x>0,'parallel pool not opened')
x = 5*x;

for i = 3:5
    a1 = rand(x:(x+i-1));
    a2 = distributed(a1);
    
    p = perms(1:i);
    for j = 1:length(p)
        d = randi(i);
        b1 = permute(a1,p(j,:));
        b2 = DistPermute(a2,p(j,:),d);
        assertEqual(b1,gather(b2));
    end
    fprintf('%i dimensions complete\n', i)
end

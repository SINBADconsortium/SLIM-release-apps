function testIterator
%TESTITERATOR - Tests basic iteration functionality of the dimTreeIterator
n = 10; d = 5; dims = n*ones(1,d);
kint = 10; kleaf = 5;
dimTree = dimensionTree(dims,kleaf,kint);

itr = dimTree.iterator('down');

assertFalse(itr.isFinished());

while itr.advance()
   
end

assertTrue(itr.isFinished());
assertTrue(~itr.advance());


end


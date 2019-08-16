function result = sinc(x)
  result = ones(size(x));
  i = (x~=0);
  if any(i(:))
    t = pi * x(i);
    result(i) = sin (t) ./ t;
  end
end

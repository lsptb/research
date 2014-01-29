function retval=file_exists(filename)
%function retval=file_exists(filename)

[status,message,messageid]=fileattrib(filename);
if (status==0)
  retval = 0;
else
  retval = 1;
end
return;

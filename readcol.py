def readcol(filename,delimiter,comment,column):
    # This function reads a text file in and returns the data in a column into a vector
    # The data is returned in the form of strings so you can cast it to whatever you need it to be
    # The delimiter of 'whitespace' will separate the lines with whitespace compressing all of the whitespace
    # into a single character. Column starts with zero.
    file = open(filename,'r')
    lines=file.readlines()
    num_lines=len(lines)
    file.close()
    data=[]
    num_col=0
    for line in lines:
        if not line[0]==comment and not line[0]=='\n':
            if delimiter=='whitespace':
                line_arr=line.split()
                num_col=max(num_col,len(line_arr))
            else:
                line_arr=line.split(delimiter)
                num_col=max(num_col,len(line_arr))
                
    for line in lines:
        if not line[0]==comment  and not line[0]=='\n':
            if delimiter=='whitespace':
                line_arr=line.split()
                
                if len(line_arr) < num_col:
                    for i in range(1,num_col-len(line_arr)):
                        line_arr.append('')
                
                data.append(line_arr[column])
            else:
                line_arr=line.split(delimiter)
             
                if len(line_arr) < num_col:
                    for i in range(1,num_col-len(line_arr)):
                        line_arr.append('')
                
                data.append(line_arr[column])
    return data

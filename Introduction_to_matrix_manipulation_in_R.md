
# Introduction to matrix manipulation in R

- A variable can hold a single value

```{r}
x <- 1; # This variable has name 'x'. It is assigned to hold the number 1
print(x);
```
- A vector holds an ordered list of values 

```
x <- 1:10; #In this case, x is a vector of cosecurtive natural numbers 1, 2, 3,..., 10
print(x);
```

- We can perform simple arithmetic operations on a vector 

```
y <- x^2 - 2*x + 4;  # Simple algebraic operations can be applied to all the values in a vector. In this case, y is also a vector of 10 values.
print(y);
plot(x, y, type="b", main="Relationship between x and y");   # The 'plot' function generates a 2-dimensional x-y plot 

```

- A matrix holds an ordered arrangement of n x m numbers (n is number of rows, and m is number of columns) 

```
x <- matrix(1:15, nrow=5, ncol=3)  #x is assigned to be a matrix consisting of 2 rows and 5 columns, containing the values 1, 2, 3,..., 15
print(x);
```  

- A matrix can be manipulated   

```
x.times.2 <- x*2; # this operation multiplies all values in x by 2
print(x.times.2);
x.plus.1 <- x+1; # this operation adds all values in x by 1
print(x.plus.1);
x.transposed <- t(x);   # this operation transposes a n x m matrix to a m x n matrix
print(x.transposed);

```  
  
  
- You can access different elements of a matrix  

```
print(x);  # Show the entire matrix
print(x[2,2]);  # Show the element in the second row and second column
print(x[1:2,2:3]); # Show the elements in row 1 and 2, and columns 2 and 3. The output is a matrix
print(x[1,]);  # Show the elements in the 1st ro. The output is a vector
print(x[,3]); # Show the elements in the 3rd column. The output is a vector
print(x[-1,]); # This operation show the elements in all rows and columns, EXCEPT the 1st row

```

- You can add row and column names to a matrix, and access elements by these names

```
rownames(x) <- paste("Row",1:5,sep="");
colnames(x) <- paste("Col", LETTERS[1:3], sep="");
print(x);
print(x["Row3",]);   
print(x[,"ColB"]);
print(x["Row3", "ColB"]);

```


- Generate summary statistics of a matrix 
```
apply(x, 1, sum); #this operation compute the sum of all values of each row (e.g., dimension 1)
apply(x, 2, mean); #this operation computer the arithmetic means of each column (e.g., dimension 2)

```


- You can use ? in front of any function to display it's usage  
```
?rownames
```




*************************
Learning objectives:
1) Understand the different data types: variable, vector and matrix.
2) Know how to enter values into and data type and recover a value from a variable, vector and matrix.
3) Perform simple arithmetic operations and manipulations in a variable, vector and matrix.
***************************




 
 
 
 
 

const arr = [1, 2, [3, 4]];
const arrCopy = [...arr];
arrCopy[2].push(5);
console.log(arr[2]); 
console.log(arrCopy[2]); 


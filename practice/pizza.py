import numpy as np
import copy
import time
import os
import sys
#%debug

#pool = Pool(processes=4)

class piz:        
    def importData(self, filename):
        """
        Header: Rows, Coloumns, Number of Ingredients, Max cells per slice then loop over "pizza", T - Tomatos, M - Mushrooms
        """
        with open(filename, 'r') as fin:
            line = fin.readline().split()
            R, C, L, H = [int(i) for i in line]
            pizza = []
            for i in range(R):
                pizza.append(np.array(list(fin.readline().strip())))
            return R, C, L, H, np.matrix(pizza)

    def export_data(self, filename, slices):
        """
        Format:
        3        3 slices.
        0 0 2 1  First slice between rows (0,2) and columns (0,1).
        0 2 2 2  Second slice between rows (0,2) and columns (2,2).
        0 3 2 4  Third slice between rows (0,2) and columns (3,4).
        [[(0, 1), (0, 0)], [(0, 4), (0, 2)], [(1, 4), (1, 2)], [(2, 4), (2, 2)], [(3, 4), (3, 0)], [(1, 1), (2, 1)], [(1, 0), (2, 0)]]
        [0, 0, 1, 2, 3, 1, 1]
        [0, 0, 1, 2, 3, 2, 2]
        [1, 4, 4, 4, 4, 1, 0]
        [0, 2, 2, 2, 0, 1, 0]

        """

        out = np.matrix([[item[0][0][0] for item in slices],
                        [item[0][1][0] for item in slices],
                        [item[0][0][1] for item in slices], 
                        [item[0][1][1] for item in slices]]).T

        with open(str(filename.split('.')[0])+".out", 'w') as f:
            f.write(str(len(slices))+"\n")
            #f.write(bytes("SP,"+lists+"\n","UTF-8"))
            #Used this line for a variable list of numbers
            np.savetxt(f, out.astype(int), fmt='%i' , delimiter=" ")
        
    def printVars(self):
        print("Rows: " + str(R))
        print("Columns: " + str(C))
        print("Minimum number of each ingredient in a slice: " + str(L))
        print("Minimum slice size: " + str(L*2))
        print("Maximum slice size: " + str(H))
        print("Maximum points: " + str(R*C))
        print(pizza)

    def oneStep(self, rows, cols, min_sl, max_sl, matrix, slices):
        #small time needed for one run 10^-5s
        #medium time needed for one run 8*10^-5s
        
        self.pizzaSlice = []
        self.tomatos = 0
        self.mushrooms = 0
        self.matrix = matrix
        self.cellsInSlice = 2*max_sl
        k = 0
        
        while self.cellsInSlice > max_sl: #max 4 runs
            
            # generate random rectangle
            row0 = np.random.randint(0, rows)
            col0 = np.random.randint(0, cols)
            
            if self.matrix[row0, col0] == 2:
                k += 1
            
            if k >= 10:
                return [], -1
            
            #print row0
            #print col0
            
            row1 = -1
            col1 = -1
            
            while not (row1 >= 0 and col1 >= 0 and row1 < rows and col1 < cols): #max 4 runs
                direction = np.random.randint(0, 4)
                if direction == 0:
                    row1 = np.random.randint(row0,row0+np.random.randint(1,max_sl))
                    col1 = np.random.randint(col0,col0+np.random.randint(1,max_sl))
                if direction == 1:
                    row1 = np.random.randint(row0,row0+np.random.randint(1,max_sl))
                    col1 = np.random.randint(col0-np.random.randint(1,max_sl),col0)
                if direction == 2:
                    row1 = np.random.randint(row0-np.random.randint(1,max_sl), row0)
                    col1 = np.random.randint(col0-np.random.randint(1,max_sl), col0)
                if direction == 3:
                    row1 = np.random.randint(row0-np.random.randint(1,max_sl),row0)
                    col1 = np.random.randint(col0,col0+np.random.randint(1,max_sl))
                #print "dir = " + str(direction)
                #print row1
                #print col1
                
            drow = [min(row0, row1), max(row0, row1)]
            dcol = [min(col0, col1), max(col0, col1)]
            
            #print dcol
            #print drow

            self.cellsInSlice = ((drow[1] - drow[0] + 1) * (dcol[1] - dcol[0] + 1))
        
        for i in range(drow[0], drow[1] + 1): #max max_sl runs
            for j in range(dcol[0], dcol[1] + 1):
                if self.matrix[i,j] == 2: # check if there is an already sliced part in the rectangle
                    return [], 0
                if self.matrix[i,j] == 0: # check how much tomatos in rectangle
                    self.tomatos += 1
                if self.matrix[i,j] == 1: # check how much mushrooms in rectangle
                    self.mushrooms += 1
                    
                    
        if self.tomatos >= min_sl and self.mushrooms >= min_sl: # check if enough mushrooms and tomatos in rect
            self.pizzaSlice.append([drow, dcol])
            for k in range(drow[0], drow[1] + 1):
                for l in range(dcol[0], dcol[1] + 1): # overwrite all cells in slice with the already sliced marker
                    self.matrix[k,l] = 2
        
            return self.pizzaSlice, self.cellsInSlice
        else:
            return [], 0

    def oneStep_XY(self, rows, cols, min_sl, max_sl, matrix, slices, x, y):
        #small time needed for one run 10^-5s
        #medium time needed for one run 8*10^-5s
        
        self.pizzaSlice = []
        self.tomatos = 0
        self.mushrooms = 0
        self.matrix = matrix
        self.cellsInSlice = 2*max_sl
        
        while self.cellsInSlice > max_sl:
            # generate random rectangle
            row0 = y
            col0 = x
            
            #print row0
            #print col0
            
            row1 = -1
            col1 = -1
            
            while not (row1 >= 0 and col1 >= 0 and row1 < rows and col1 < cols):
                direction = np.random.randint(0, 4)
                if direction == 0:
                    row1 = np.random.randint(row0,row0+np.random.randint(1,max_sl))
                    col1 = np.random.randint(col0,col0+np.random.randint(1,max_sl))
                if direction == 1:
                    row1 = np.random.randint(row0,row0+np.random.randint(1,max_sl))
                    col1 = np.random.randint(col0-np.random.randint(1,max_sl),col0)
                if direction == 2:
                    row1 = np.random.randint(row0-np.random.randint(1,max_sl), row0)
                    col1 = np.random.randint(col0-np.random.randint(1,max_sl), col0)
                if direction == 3:
                    row1 = np.random.randint(row0-np.random.randint(1,max_sl),row0)
                    col1 = np.random.randint(col0,col0+np.random.randint(1,max_sl))
                #print "dir = " + str(direction)
                #print row1
                #print col1

            drow = [min(row0, row1), max(row0, row1)]
            dcol = [min(col0, col1), max(col0, col1)]
            
            if not dcol[1]+1 >= cols:
                x = dcol[1]+1
            else:
                x = 0
                if not drow[1]+1 >= rows:
                    y = drow[1]+1
                else:
                    y = 0
            
            #print dcol
            #print drow

            self.cellsInSlice = ((drow[1] - drow[0] + 1) * (dcol[1] - dcol[0] + 1))
        
        #if self.cellsInSlice <= max_sl: # check if rectangle has less than allowed cells
        for i in range(drow[0],drow[1] + 1):
            for j in range(dcol[0],dcol[1] + 1):
                if self.matrix[i,j] == 2: # check if there is an already sliced part in the rectangle
                    return [], 0
                if self.matrix[i,j] == 0: # check how much tomatos in rectangle
                    self.tomatos += 1
                if self.matrix[i,j] == 1: # check how much mushrooms in rectangle
                    self.mushrooms += 1
            
        if self.tomatos >= min_sl and self.mushrooms >= min_sl: # check if enough mushrooms and tomatos in rect
            self.pizzaSlice.append([drow, dcol])
            for k in range(drow[0], drow[1] + 1):
                for l in range(dcol[0], dcol[1] + 1): # overwrite all cells in slice with the already sliced marker
                    self.matrix[k,l] = 2
        
            return self.pizzaSlice, self.cellsInSlice
        else:
            return [], 0

    def resetPizza(self, inpizza):
        return copy.deepcopy(inpizza)
    
    def main(self, pizza, R, C, L, H, prec):
        try:
            start = time.time()
            self.pizzaSlice = []
            self.slices = []
            self.sliceing = 0
            maxSlices = 0
            self.pizza = pizza
            self.pizzatmp = piz.resetPizza(self.pizza)

            #for p in range(0,10):
            k = 0

            if prec <= 1:
                maxSlices = 0
                while self.sliceing < int(R*C*prec):
                    k += 1
                    self.pizzaSlice = []
                    self.slices = []
                    self.sliceing = 0
                    self.pizzatmp = piz.resetPizza(self.pizza)
                    for i in range(0,R*C*8):

                        self.pizzaSlice, self.cellNr = piz.oneStep(R, C, L, H, self.pizzatmp, self.slices)
                        
                        if self.cellNr == -1:
                            break
                        if self.pizzaSlice != []:
                            self.slices.append(self.pizzaSlice)
                            self.sliceing += self.cellNr
                        if self.sliceing >= R*C:
                            break

                    if self.sliceing > maxSlices:
                        maxSlices = self.sliceing
                        print maxSlices,
                        self.slices = self.slices
            else:
                maxSlices = 0
                for p in range(0,prec):

                    self.pizzaSlice = []
                    self.slices = []
                    self.sliceing = 0
                    self.pizzatmp = piz.resetPizza(self.pizza)

                    for i in range(0,R*C*8):
                        k += 1

                        #start = time.time()
                        self.pizzaSlice, self.cellNr = piz.oneStep(R, C, L, H, self.pizzatmp, self.slices)
                        #print "onestep time:" + str(time.time() - start)

                        if self.cellNr == -1:
                            break
                        if self.pizzaSlice != []:
                            self.slices.append(self.pizzaSlice)
                            self.sliceing += self.cellNr
                        if self.sliceing >= R*C:
                            break

                    if self.sliceing > maxSlices:
                        maxSlices = self.sliceing
                        print "Max slices found " + str(maxSlices) + " at iteration " + str(k) + " of " + str(R*C*R/10*prec)
                        self.slices = self.slices

        except KeyboardInterrupt:
            return maxSlices, self.slices
            try:
                sys.exit(0)
            except SystemExit:
                os._exit(0)
        
        end = time.time()
        print "\n"
        print(end - start)
        print(maxSlices)
        #print(slices)
        return maxSlices, self.slices

filename = "big.in"
lastMax = 0
#np.random.seed(seed=5)

piz = piz()
R, C, L, H, pizza = piz.importData(filename)

piz.printVars()

pizzaDual = piz.resetPizza(pizza)
pizzaDual[pizzaDual == 'T'] = 0
pizzaDual[pizzaDual == 'M'] = 1
pizzaDual = np.matrix(pizzaDual.astype(int))
print pizzaDual

n, slices = piz.main(pizzaDual, R, C, L, H, 2)
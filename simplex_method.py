# python3
from sys import stdin
import itertools
EPS = 1e-4

class Position:
    def __init__(self, column, row):
        self.column = column
        self.row = row

# Takes an int "n" that defines how many total inequalty equations there are
# Takes an int "m" that defines how many total variables there are
# Takes a list of lists "a" of size n x m that contains the inequalities.
# Takes a list "b" of size n of the maximums of each inequality
# Takes a list "c" that represents the optimization function to maximize
# For example:
# x + y - 3z <= 10
# -5x + 10y <= 50
# 3x - 2y -4z <= 9
# Maximize: x + 6y -3z
# n, m = 3, 3
# a = [[1,1,-3],[-5,10,0],[3,-2,-4]]
# b = [10, 50, 9]
# c = [1, 6, -3]
def ReadEquation():
    n, m = map(int, input().split())
    a = []
    for row in range(n):
        a.append(list(map(float, input().split())))
    b = list(map(float, input().split()))
    c = list(map(float, input().split()))
    return a, b, c, n, m

#A tableau is created of size (m + n + 1) x (n + 1). Each inequality is placed in a row in the tableau and a slack variable is added to
#each inequality. A slack variable in inequality equations are added to transform the equations to an equalities.
#In the provided example, the tableau would look as such for the standard format (non 2 phase structure):
# [ 1,  1, -3, 1, 0, 0, 10]
# [-5, 10,  0, 0, 1, 0, 50]
# [ 3, -2, -4, 0, 0, 1,  9]
# [-1, -6,  3, 0, 0, 0,  0]
def CreateTableau(a, b, c, n, phase_one_optimization):
  tableau = []
  phase_one_row = [0] * (len(c) + n + 2)
  for i in range (n):
    #For phase 1 optimization in a the 2-phase simplex method, any inequalities that have a solution less than zero will be flipped
    #For example 2x -3y <= -10 --> -2x + 3y <= 10
    #The slack variable is added to the tableau with value -1 to account for the flip in sign
    #Recall that the two phase approach will ONLY occur if an optimal solution was not initially found
    if phase_one_optimization and b[i] < 0:
      slack_variables = [0] * n
      slack_variables[i] = -1.0
      tableau_row = [-1*x for x in a[i]] + slack_variables + [-1 * b[i]]
      tableau.append(tableau_row)
      phase_one_row = [a + b for a, b in zip(phase_one_row, tableau_row)]
    else:
      slack_variables = [0] * n
      slack_variables[i] = 1.0
      tableau_row = a[i] + slack_variables + [b[i]]
      tableau.append(tableau_row)
  final_row = [-1*x for x in c] + [0] * n + [0]
  tableau.append(final_row)
  return tableau, phase_one_row

#Bland's Rule will be used for selecting the pivot element:
# 1. Choose the leftmost column that is negative
# 2. Among the rows, choose the one with the lowest ratio between the right-hand side of the tableau (value b) and the column coefficint where the coefficient is greater than zero.
# 2 (cont). If the minimum ratio is shared by several rows, choose the row with the lowest column variable (basic variable) in it.
# For #2, the algorithm doesn't just take the min low ratio because the special case of multiple minimums must be taken into account. Additionally, bland's rule calls for the lowest-
# numbered basic variable which is different than the lowest index. For this, the algorithm keeps account of the basic variable in the "slack_rows" list
def SelectPivotElement(a, m, slack_rows, phase_one_optimization, phase_one_row):
    pivot_element = Position(0, 0)
    no_solution = False
    if phase_one_optimization:
      pivot_element.column = phase_one_row.index(max(phase_one_row[:-1]))
    else:
      pivot_element.column = a[len(a)-1][:-1].index(min(a[len(a)-1][:-1])) #Choose minimum based on first negative smallest index
    ratios = []
    if pivot_element.column != None:
      for r in range(len(a)-1):
        if a[r][pivot_element.column] > 0:
            ratios.append(abs(a[r][-1] / a[r][pivot_element.column]))
        else:
          ratios.append(float("inf"))
      if all(i == float("inf") for i in ratios):
        no_solution = True
      row_min = min(ratios)
      row_min_indicies = [i for i,x in enumerate(ratios) if x == row_min]
      #take into account the case of equal minimums in rows. According to Bland's rule, choose least variable
      if (len(row_min_indicies) > 1):
        least_variable = []
        for j in row_min_indicies:
          least_variable.append(slack_rows[j])
        pivot_element.row = slack_rows.index(min(least_variable))
      else:
        pivot_element.row = row_min_indicies[0]
    else:
      no_solution = True
    return no_solution, pivot_element

#Process Pivot has been optimized to insert the value 0 in the pivot element column. With real number manipulation, this helps ensure within the tableau the reduction to zero is always
#clear and not approximated with the episilon functions.
def ProcessPivotElement(a,pivot_element, phase_one_optimization, phase_one_row):
    pri_mult = a[pivot_element.row][pivot_element.column] #primary multiplier from pivot element
    a[pivot_element.row] = [n / pri_mult for n in a[pivot_element.row]] #make primary element have a value of 1
    a[pivot_element.row][pivot_element.column] = 1.0
    for i in range(len(a)):
      if i != pivot_element.row:
        sec_mult = a[i][pivot_element.column] #secondary multiplier from row being updated
        pri_row = [j * sec_mult for j in a[pivot_element.row]]
        a[i]= [a- b for a, b in zip(a[i], pri_row)]
        a[i][pivot_element.column] = 0
    if phase_one_optimization:
      sec_mult = phase_one_row[pivot_element.column] #secondary multiplier from row being updated
      pri_row = [j * sec_mult for j in a[pivot_element.row]]
      phase_one_row = [a- b for a, b in zip(phase_one_row, pri_row)]
      phase_one_row[pivot_element.column] = 0
    return a, phase_one_row

#Solves a linear program inequality use a tableau through the simplex method
#The algorithm will first attempt to solve the tableau assuming a basic feasible solution has been provided. If the tableau provided by the first attempt leads to an invalid solution,
#meaning one of the inequality equations is violated with one of the values in the initial optimal values, then the algorithm will create a new tableau 
#and proceed to a two-phase simplex method approach.
def SolveEquation(a, b, c, n, m):
    if all(i <= 0 for i in c) and all(i >= 0 for i in b):
      return [0] * m
    tableau, phase_one_row = CreateTableau(a, b, c, n, False)
    ans, phase_one_answer = solveTableau(tableau, a, b, m, n, False, phase_one_row)
    #break immediately if the tableau reduced to 
    if ans == [-1] or ans == [float("inf")]:
      return ans
    invalid_answer = valid_answer(ans, a, b, m, n)
    #Proceed to a two-phase simplex approach if one of the variables in the optimial solution violates an inequality equation
    if invalid_answer:
      tableau, phase_one_row = CreateTableau(a, b, c, n, True)
      ans, phase_one_answer = solveTableau(tableau, a, b, m, n, True, phase_one_row)
      phase_one_answer_invalid = valid_answer(phase_one_answer, a, b, m, n)
      if ans == [-1] or ans == [float("inf")]:
        return ans
      invalid_answer = valid_answer(ans, a, b, m, n)
    if invalid_answer:
      if not phase_one_answer_invalid:
        return phase_one_answer
      else:
        return [-1]
    return ans

def valid_answer(ans, a, b, m, n):
  invalid_answer = False
  for i in range(n):
        valid_ans = 0
        for j in range(m):
          valid_ans += a[i][j] * ans[j]
        if epsilon_greater_than(valid_ans, b[i]):
          invalid_answer = True
  if not all(epsilon_greater_than_equal_to(i, 0) for i in ans):
    invalid_answer = True
  return invalid_answer

def solveTableau(tableau, a, b, m, n, phase_one_optimization, phase_one_row):
  slack_rows = list(range(m,n+m))
  phase_one_complete = False
  phase_one_answer = [0] * m
  while (phase_one_optimization or not all(epsilon_greater_than_equal_to(i, 0) for i in tableau[len(tableau)-1][:-1])):
    if phase_one_optimization and all(epsilon_less_than_equal_to(k, 0) for k in phase_one_row[:-1]):
      phase_one_optimization = False
      phase_one_complete = True
      phase_one_answer = determine_answer(tableau, slack_rows)
      if all(epsilon_greater_than_equal_to(i, 0) for i in tableau[len(tableau)-1][:-1]):
        break
    no_solution, pivot_element = SelectPivotElement(tableau, m, slack_rows, phase_one_optimization, phase_one_row)
    if no_solution:
      if phase_one_complete:
        return [-1], phase_one_answer
      else:
         return [float("inf")], phase_one_answer
    slack_rows[pivot_element.row] = pivot_element.column
    tableau, phase_one_row = ProcessPivotElement(tableau, pivot_element, phase_one_optimization, phase_one_row)
  return determine_answer(tableau, slack_rows), phase_one_answer

def determine_answer(tableau, slack):
  ans = [0] * m
  for i in range(n+m):
    if i < m and i in slack:
      index = slack.index(i)
      ans[i] = tableau[index][-1]
    elif i not in slack and tableau[-1][i] == 0:
      for j in range(n-1):
        if tableau[j][i] > 0:
          return [-1]
    elif i < m:
      ans[i] = 0
  return ans

#Measure equality or inequality of two real numbers through an epsilon value EPS
def epsilon_greater_than(a, b):
  return ((a > b) and not isclose(a, b))

def epsilon_greater_than_equal_to(a, b):
  return ((a > b) or isclose(a, b))

def epsilon_less_than(a, b):
  return ((a < b) and not isclose(a, b))

def epsilon_less_than_equal_to(a, b):
  return ((a < b) or isclose(a, b))

def isclose(a, b):
    return abs(a-b) <=EPS

def PrintColumn(column):
    size = len(column)
    if size == 1 and column[0] == -1:
      print("No solution")
    elif size == 1 and column[0] == float("inf"):
      print("Infinity")
    else:
      print("Bounded solution")
      print(' '.join(list(map(lambda x : '%.18f' % x, column))))

if __name__ == "__main__":
    a, b, c, n, m = ReadEquation()
    solution = SolveEquation(a, b, c, n, m)
    PrintColumn(solution)
    exit(0)
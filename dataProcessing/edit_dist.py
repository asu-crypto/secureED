
def fetchComparisions(a, b):
    E = [[1 for _ in range(len(b))] for _ in range(len(a))]
    for i in range(len(a)):
        for j in range(len(b)):
            if a[i] == b[j]:
                E[i][j] = 0
    return E

def edit_distance_R(a, b):

    m, n = len(a), len(b)
    E = fetchComparisions(a, b)

    R = [i for i in range(n)]

    for i in range(m):
        temp = R[0]
        R[0] = min(i+E[i][0], R[0]+1)
        for j in range(1,n):
            temp2 = R[j]
            R[j] = min(temp+E[i][j], R[j-1]+1, R[j]+1)
            temp = temp2

    return R[-1], R

def edit_distance_D(a, b):

    m, n = len(a), len(b)
    E = fetchComparisions(a, b)

    D = [[0 for _ in range(n+1)] for _ in range(m+1)]
    for i in range(m+1):
        D[i][0] = i
    for j in range(1, n+1):
        D[0][j] = j

    for i in range(m):
        for j in range(n):
            D[i+1][j+1] = min(D[i][j]+E[i][j], D[i][j+1]+1, D[i+1][j]+1)
    return D[-1][-1], D
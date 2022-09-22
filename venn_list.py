
#return_type=tuple(([1,2],[2,3],[3,4]))

def times_two(x=5):
    return x*2

times_two() # = ERROR
times_two(3) # = 6

def venn_list(return_type=None,*args):
    if return_type is None:
        return_type = tuple(([1,2],[2,3],[3,4]))
    print(return_type,type(return_type))
    print(args)
    print(list(return_type)+list(args))
    args=list(return_type)+list(args)
    for arg in args:
        print(arg)
        if not isinstance(arg,list):
            raise TypeError("Only Lists are allowed")
    intersection_l=set.intersection(*map(set,args))
    union_l=set.union(*map(set,args))
    return tuple(intersection_l),tuple(union_l)

l0=[[1,3],[2,3],[3,4]]
i0,u0=venn_list(None,*l0)
print(i0,u0)

l=[[1,2,9],[2,3,4],[3,4,5]]
i,u=venn_list(None,*l)
print(i)
print(u)
print(list(map(type,[i,u])))

print("Error:")
l1=["A",[1,2,3]]
#l1=[[1,2,3],[1,2,3]]
i1,u1=venn_list(None,*l1)

def venn_list1(return_type=([1,2],[2,3],[3,4],),*args):
    print(return_type)
    print(args)

l2=[[1,2,9],[2,3,4],[3,4,5]]
i2,u2=venn_list1(*l2)
print(i2)
print(u2)

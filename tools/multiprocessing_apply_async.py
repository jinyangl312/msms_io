import multiprocessing


def calculate(func, args):
    return func(*args)


def calculatestar(args):
    return calculate(*args)


def apply_async_with_batch(processer_num, task_generator, callback, batch_size):
    try:
        while True:
            pool = multiprocessing.Pool(processer_num)
            for i in range(batch_size):
                pool.apply_async(calculate, next(task_generator),
                                 callback=callback)  # unordered
                # pool.apply(calculate, next(task_generator))  # unordered
            pool.close()
            pool.join()
    except Exception as e:
        # print(repr(e))
        if 'pool' in vars() and isinstance(pool, multiprocessing.pool.Pool):
            pool.close()
            pool.join()

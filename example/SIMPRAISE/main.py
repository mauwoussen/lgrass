# coding: utf8
# '''
# Created on 21/04/2020
#
# @author: modelisation - TR
# '''
import multiprocessing
import time
import pandas as pd
import lgrass_batch_simpraise as batch


if __name__ == '__main__':
    timing = time.time()
    plan = pd.read_csv("example/SIMPRAISE/inputs/plan_simulation.csv", sep=',')
    genet_path = "example/SIMPRAISE/modelgenet"

    for i in range(1,4):
        batch.runlsystem(plan, id_scenario=i, genet_path=genet_path, id_gener=1)

    ###    Utilisation de Pool multiprocess    ###

    # multiprocessing.freeze_support()
    # CPUnb = 3  # int(multiprocessing.cpu_count())  # nombre de processeurs utilises
    # pool = multiprocessing.Pool(processes=CPUnb)
    # for j in range(1, 4):
    #     pool.apply_async(batch.runlsystem, args=(plan, j, 1))
    # pool.close()
    # pool.join()
    # print('Global execution time : ', time.time() - timing)

    ###   Utilisation de Process multiprocess   ###

    # q = multiprocessing.Queue()
    # p = multiprocessing.Process(target=batch.runlsystem, args=(plan, 16, 1, q))
    # p.start()
    # p.join()
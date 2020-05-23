for nCores=1 % 6:-1:1
    nCores
    maxNumCompThreads(nCores)
    compare_all
    save(sprintf('result_cores%d',nCores),'nCores','timeTable','timeTables');
end

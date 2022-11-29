from aPhyloGeo import geneticPipeline, climaticPipeline

def main ():
    climaticTrees = climaticPipeline()
    geneticPipeline(climaticTrees)

if __name__ == "__main__":
    main()
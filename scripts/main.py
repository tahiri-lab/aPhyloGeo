from aPhyloGeo import geneticPipeline, climaticPipeline

if __name__ == "__main__":
    climaticTrees = climaticPipeline()
    geneticPipeline(climaticTrees)
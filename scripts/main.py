from aPhyloGeo import geneticPipeline, climaticPipeline
titleCard=r"""                               
        ____    __               ___                                           
       /\  _`\ /\ \             /\_ \                                          
   __  \ \ \L\ \ \ \___   __  __\//\ \     ___      __      __    ___          
 /'__`\ \ \ ,__/\ \  _ `\/\ \/\ \ \ \ \   / __`\  /'_ `\  /'__`\ / __`\        
/\ \L\.\_\ \ \/  \ \ \ \ \ \ \_\ \ \_\ \_/\ \L\ \/\ \L\ \/\  __//\ \L\ \       
\ \__/.\_\\ \_\   \ \_\ \_\/`____ \/\____\ \____/\ \____ \ \____\ \____/       
 \/__/\/_/ \/_/    \/_/\/_/`/___/> \/____/\/___/  \/___L\ \/____/\/___/        
                              /\___/                /\____/                    
                              \/__/                 \_/__/                     
"""

if __name__ == "__main__":
    print(titleCard+"\n")
    climaticTrees = climaticPipeline()
    geneticPipeline(climaticTrees)
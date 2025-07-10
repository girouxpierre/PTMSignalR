from neo4j import GraphDatabase

def requestPhospho(uri, user , password, queryCataSolo, queryCataComp, queryGenePin, queryGenePout, outputPrefix):
    """
    launch request on the neo4j reactome db to get genes which controls phospho

    """
    
    driver = GraphDatabase.driver(uri, auth=(user, password))
    
    with driver.session() as session:
        result1 = session.run(queryCataSolo)
        pathName1 = outputPrefix+'_CataSolo.tsv'
        with open(pathName1, 'w', encoding='utf-8') as out1:
            out1.write("reaction\tcatalysers\n")
            
            for record in result1:
                out1.write(f"{record['reaction']}\t{record['catalysers']}\n")

        result2 = session.run(queryCataComp)
        pathName2 = outputPrefix+'_CataComp.tsv'
        with open(pathName2, 'w', encoding='utf-8') as out2:
            out2.write("reaction\tcatalysers\n")
            
            for record in result2:
                out2.write(f"{record['reaction']}\t{record['catalysers']}\n")   

        result3 = session.run(queryGenePin)
        pathName3 = outputPrefix+'_GenePin.tsv'
        with open(pathName3, 'w', encoding='utf-8') as out3:
            out3.write("reaction\tgeneName\tmodifiedPositions\n")
            
            for record in result3:
                out3.write(f"{record['reaction']}\t{record['geneName']}\t{record['modifiedPositions']}\n")

        result4 = session.run(queryGenePout)
        pathName4 = outputPrefix+'_GenePout.tsv'
        with open(pathName4, 'w', encoding='utf-8') as out4:
            out4.write("reaction\tgeneName\tmodifiedPositions\n")
            
            for record in result4:
                out4.write(f"{record['reaction']}\t{record['geneName']}\t{record['modifiedPositions']}\n") 

    driver.close()         


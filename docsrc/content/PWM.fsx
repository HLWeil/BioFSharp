(*** hide ***)
// This block of code is omitted in the generated HTML documentation. Use 
// it to define helpers that you do not want to show in the documentation.
#I @"../../bin/BioFSharp/net47/"
#I @"../../bin/BioFSharp.BioDB/net47/"
#I @"../../bin/BioFSharp.ImgP/net47"
#I @"../../bin/BioFSharp.IO/net47/"
#I @"../../bin/BioFSharp.Parallel/net47/"
#I @"../../bin/BioFSharp.Stats/net47/"
#I @"../../bin/BioFSharp.Vis/net47/"


#r @"../../packages/formatting/FSharp.Plotly/lib/netstandard2.0/FSharp.Plotly.dll"

#r "netstandard.dll"
#r "BioFSharp.dll"
#r "BioFSharp.IO.dll"
#r "FSharpAux.dll"
#r "FSharpAux.IO.dll"
#r @"C:\Users\weil\source\repos\FSharp.Stats\src\FSharp.Stats\bin\Release\netstandard2.0\FSharp.Stats.dll"


//#r "FSharpAux.IO.dll"

(**
<table class="HeadAPI">
<td class="Head"><h1>Position Specific Scoring Matrices</h1></td>
<td class="API">
    <a id="APILink" href="https://csbiology.github.io/BioFSharp/reference/biofsharp-pwm.html" >&#128194;View module documentation</a>
</td>

**)

open BioFSharp
open BioFSharp.IO
open Nucleotides
open PositionSpecificScoringMatrices
open FSharp.Plotly
open FSharp.Stats
open FSharp.Stats.ML
open FSharp.Stats.ML.Unsupervised

//fsi.AddPrinter(BioFSharp.IO.FSIPrinters.prettyPrintPFM)
//fsi.AddPrinter(BioFSharp.IO.FSIPrinters.prettyPrintPWM)
//fsi.AddPrinter(BioFSharp.IO.FSIPrinters.prettyPrintPPM)

fsi.AddPrinter(BioFSharp.IO.FSIPrinters.prettyPrintBioItemFrequencies)

let sequencesForMotif1 = 
    [|
        "GTTCAG"
        "ATTCAG"
        "ATCTAG"
        "GTTAAG"
        "ATTAAG"   
    |]
    |> Array.map BioArray.ofNucleotideString

let sequencesForMotif2 = 
    [|
        "GTTCAGA"
        "AGTCAGA"
        "GTCTAAT"
        "GTTAAGA"
        "ATTAAGG"   
    |]
    |> Array.map BioArray.ofNucleotideString

let sequencesForMotif3 = 
    [|
        "TTCATCC"
        "TGCATCC"
        "TGCATTC"
        "AGTATCC"
        "ATTATCC"   
    |]
    |> Array.map BioArray.ofNucleotideString 


let nucleotides = [|Nucleotide.A;Nucleotide.G;Nucleotide.C;Nucleotide.T|] 

let backGround = 
    BackgroundFrequencies.createFCVOfSequence nucleotides
    |> BackgroundFrequencies.createPCVOfFCV nucleotides  0.



let createPWMOfSequences sequences =
    sequences
    |> Array.map (fun s -> 
        PositionMatrix.createPFMOfSingleSequence s
    )
    |> PositionMatrix.fusePositionFrequencyMatrices sequences.[0].Length
    |> PositionMatrix.createPPMOfPFM 5 nucleotides 0.
    |> PositionMatrix.createPositionWeightMatrix nucleotides backGround

let motif1 = createPWMOfSequences sequencesForMotif1 
let motif2 = createPWMOfSequences sequencesForMotif2
let motif3 = createPWMOfSequences sequencesForMotif3


open HierarchicalClustering


/// Creates a phylogenetic tree from the hierarchical cluster. The nameF is used to name branches based on their id.  
let toPhylTree (nameF : int -> 'T) (cluster: Cluster<'T>) =
    let rec loop upperDistance cluster =
        match cluster with
        | Node (id,dist,num,c1,c2) ->
            PhylTree.Branch ((id,upperDistance, nameF id), [loop dist c1; loop dist c2])
        | Leaf (id,x,name) -> PhylTree.Branch ((id,upperDistance,name),List.empty)
    loop 0. cluster

/// Creates a hierarchical cluster from the phylogenetic tree.
let ofPhylTree (tree: PhylTree.Node<int*float*'T>) =
    let rec loop tree =
        match tree with
        | PhylTree.Branch ((id,upperDist,name),[]) ->
            Leaf(id,1,name),1,upperDist
        | PhylTree.Branch ((id,upperDist,name),[b1;b2]) -> 
            let c1,num1,dist1 = loop b1
            let c2,num2,dist2 = loop b2
            let dist = (dist1 + dist2) / 2.
            let num = num1 + num2
            Node (id,dist,num,c1,c2),num,upperDist
        | PhylTree.Branch ((id,dist,name),l) -> 
            failwithf "Cannot transform branch of phylogenetic tree with id %i to hierarchical cluster, as it needs to have either no or 2 subbranches but here has %i" id l.Length
    let c,_,_ = loop tree
    c


Similarity.compareTwoPositionMatrices Correlation.Seq.pearson motif1 motif2

open FSharpAux

let motif4 = 
    FSharpAux.IO.FileIO.readFile "C:\Users\weil\Downloads\JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar\MA0002.2.jaspar"
    |> Seq.toList



let alphabetF = BioItemsConverter.OptionConverter.charToOptionStandardNucleotid >> Option.get

let inline valueF x = float x


//[
//PositionMatrix.createTaggedPSSM "dudeson" motif1
//PositionMatrix.createTaggedPSSM "duderony" motif2
//]
//|> BioFSharp.IO.PositionSpecificScoringMatrices.Write.writePSSMs 
//    BioFSharp.IO.PositionSpecificScoringMatrices.Write.writeJaspar id string 
//    @"E:\Lukas\testMotifsJaspar.tf"

let xs = BioFSharp.IO.PositionSpecificScoringMatrices.Read.readMultipleJasparMatrices alphabetF valueF @"E:\Lukas\testMotifsJaspar.tf" |> Seq.toArray


let readM name path =
    let ts = BioFSharp.IO.PositionSpecificScoringMatrices.Read.readCisBPMatrix alphabetF valueF name path
    ts.Matrix 
    |> PositionMatrix.positionProbabilityMatrix
    |> PositionMatrix.createTaggedPSSM ts.Tag

let motifPaths = 
    FSharpAux.IO.FileIO.filesInDir (System.IO.DirectoryInfo(@"E:\Lukas\Projekte\ArabidopsisGRN\material\Arabidopsis_thaliana_2019_07_22_1_31_pm\pwms_all_motifs"))
    |> Array.map (fun file -> file.Name.Split '.' |> Array.head, file.FullName)

let matrices =
    motifPaths
    |> Array.map (fun (name,path) -> readM name path)
    |> Array.filter (fun tm -> tm.Matrix.Length > 4)

let pearsonFixed (xs :float []) (ys :float []) = 
    if 
        Array.fold (fun (b,v1) v2 -> v1 = v2 && b,v2) (true,xs.[0]) xs |> fst
        ||
        Array.fold (fun (b,v1) v2 -> v1 = v2 && b,v2) (true,ys.[0]) ys |> fst
        then 0.
    else Correlation.Seq.pearson xs ys

Array.fold (fun (b,v1) v2 -> v1 = v2 && b,v2) (true,1.) [|2.;1.;1.|]

#time
let similarities = Similarity.comparePositionMatrices pearsonFixed matrices

let distanceFunction = 
    let similarityMap =
        similarities
        |> Array.collect (fun ((a,b),score) -> [|((a,b),score);((b,a),score)|])
        |> Map.ofArray
    fun a b -> 1. - similarityMap.[a,b]


let similarityCluster = 
    matrices
    |> Array.map (fun tm -> tm.Tag)
    |> HierarchicalClustering.generate distanceFunction HierarchicalClustering.Linker.wardLwLinker

let phylTree = 
    similarityCluster
    |> toPhylTree (fun _ -> "Node")

phylTree
|> IO.Newick.toFile (fun (id,dist,name) -> sprintf "%i_%s" id name, sprintf "%f" dist) @"E:\Lukas\Projekte\ArabidopsisGRN\material\MotifClustering.newick"

phylTree
|> PhylTree.tryGetNodeBy (fun n -> 
    match n with 
    | PhylTree.Branch ((id,dist,name),_) -> 
        id = 2977
)

HierarchicalClustering.getDistancesOfCluster similarityCluster
|> List.countBy (fun x -> x.Equals(nan))



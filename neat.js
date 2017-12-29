function Neat(population, inputSize, outputSize) {
	inputSize++; //Bias
	
	var deltaValues = {
		deltaDisjoint: 2.0,
		deltaWeights: 0.4,
		deltaThreshold: 1.0
	};
	
	var defaultRates = {
		connections: 0.25,
		link: 2.0,
		bias: 0.4,
		node: 0.5,
		enable: 0.2,
		disable: 0.4,
		step: 0.1
	};
	
	var chances = {
		mutateConnectionsChance: 0.25,
		perturbChance: 0.90,
		crossoverChance: 0.75,
		linkMutationChance: 2.0,
		nodeMutationChance: 0.50,
		biasMutationChance: 0.40,
		disableMutationChance: 0.4,
		enableMutationChance: 0.2
	}
	
	var maxNodes = 1000000;
	var staleSpecies = 15;
	//var population = p;

	var pool = {};
	
	var sigmoid = function(x) {
		return 2/(1+Math.exp(-4.9*x))-1;
	}
	
	var newInnovation = function(){
		return ++pool.innovation;
	}
	
	var newPool = function(){
		return {
			species: [],
			generation: 0,
			innovation: outputSize,
			maxFitness: 0
		};
	}
	
	var newSpecie = function(){
		return {
			topFitness: 0,
			staleness: 0,
			genomes: [],
			averageFitness: 0
		}
	}
	
	var newGenome = function(){
		return {
			genes: [],
			fitness: 0,
			adjustedFitness: 0,
			network: {},
			maxneuron: 0,
			globalRank: 0,
			mutationRates: defaultRates
		}
	}

	var copyGenome = function(genome){
		var cGenome = {
			genes: [],
			maxneuron: genome.maxneuron,
			fitness: 0,
			mutationRates: {}
		};
		for(var i = 0; i < genome.genes.length; i++){
			cGenome.genes.push(copyGene(genome.genes[i]));
		}
		
		for(var n in genome.mutationRates){
			cGenome.mutationRates[n] = genome.mutationRates[n];
		}
		
		return cGenome;
	}
	
	var basicGenome = function(){
		var genome = newGenome();
		innovation = 1;
		
		genome.maxneuron = inputSize;
		mutate(genome);
		return genome;
	}
	
	var newGene = function(){
		return {
			into: 0,
			out: 0,
			weight: 0.0,
			enabled: true,
			innovation: 0
		};
	}
	
	var copyGene = function(gene){
		var nGene = newGene();
		for(var k in gene){
			nGene[k] = gene[k];
		}
		return nGene;
	}
	
	var newNeuron = function(){
		return {
			incoming: [],
			value: 0.0
		}
	}
	
	var evaluateNetwork = function(network, inputs) {
		inputs.push(1);
		if(inputs.length != inputSize){
			console.log("Incorrect inputs	");
			return {};
		}
		
		for(var k in inputs){
			network.neurons[k].value = inputs[k];
		}
		
		var sum = 0;
		for(var k in network.neurons){
			var neuron = network.neurons[k];
			sum = 0;
			for(var kI in neuron.incoming){
				var incoming = neuron.incoming[kI];
				//console.log(incoming,network);
				var other = network.neurons[incoming.into];
				sum += incoming.weight * other.value;
			}
			
			if(neuron.incoming.length > 0){
				neuron.value = sigmoid(sum);
			}
		}
		var outputs = {};
		var bigger = {i: 0, v: 0};
		for(var i = 0; i < outputSize; i++) {
			outputs[i] = 0;
			
			if(network.neurons[maxNodes + i].value > bigger.v) {
				bigger.i = i;
				bigger.v = network.neurons[maxNodes + i].value;
			}
		}
		if(bigger.v > 0){
			outputs[bigger.i] = 1;
		}
		return outputs;
	}
	
	var crossover = function(g1,g2){
		if(g2.fitness > g1.fitness) {
			var tmp = g1;
			g1 = g2;
			g2 = tmp;
		}
		
		var child = newGenome();
		
		var innovations2 = {}
		for(var i = 0; i < g2.genes.length; i++){
			var gene = g2.genes[i];
			innovations2[gene.innovation] = gene;
		}
		
		for(var i = 0; i < g1.genes.length; i++) {
			var gene1 = g1.genes[i];
			var gene2 = innovations2[gene1.innovation];
			if(gene2 != null && Math.random(2) == 1 && gene2.enabled){
				child.genes.push(copyGene(gene2));
			} else {
				child.genes.push(copyGene(gene1));
			}
		}
		child.maxneuron = Math.max(g1.maxneuron,g2.maxneuron);
		
		for(var k in g1.mutationRates){
			child.mutationRates[k] = g1.mutationRates[k];
		}
		
		return child;
	}
	
	var randomNeuron = function(genes,nonInput) {
		var neurons = {};
		if(!nonInput){
			for(var i = 0; i < inputSize; i++){
				neurons[i] = true;
			}
		}
		
		for(var i = 0; i < outputSize; i++){
			neurons[maxNodes+i] = true;
		}
		
		for(var i = 0; i < genes.length; i++){
			if(!nonInput || genes[i].into > inputSize) {
				neurons[genes[i].into] = true;
			}
			
			if(!nonInput || genes[i].out > inputSize) {
				neurons[genes[i].out] = true;
			}
		}
		var count = 0;
		for(var k in neurons) {
			count++;
		}
		var n = parseInt(Math.random()*count);
		for(var k in neurons){
			n--;
			if(n == 0) {
				return k;
			}
		}
		return 0;
	}

	var containsLink = function(genes,link){
		for(var i = 0; i < genes.length; i++){
			var gene = genes[i];
			if(gene.into == link.into && gene.out == link.out){
				return true;
			}
		}
		return false;
	}
	
	var pointMutate = function(genome){
		var step = genome.mutationRates['step'];
		
		for(var i = 0; i < genome.genes.length; i++) {
			var gene = genome.genes[i];
			if(Math.random() < chances.perturbChance){
				gene.weight += Math.random() * step * 2 - step;
			} else {
				gene.weight = Math.random() * 4-2;
			}
		}
	}
	
	var linkMutate = function(genome, forceBias) {
		var neuron1 = randomNeuron(genome.genes,false);
		var neuron2 = randomNeuron(genome.genes,true);
		
		var newLink = newGene();
		if(neuron1 <= inputSize && neuron2 <= inputSize){
			return;
		}
		
		if(neuron2 <= inputSize) {
			var tmp = neuron1;
			neuron1 = neuron2;
			neuron2 = tmp;
		}
		
		newLink.into = neuron1;
		newLink.out = neuron2;
		
		if(forceBias){
			newLink.into = inputSize;
		}
		if(containsLink(genome.genes,newLink)) {
			return;
		}
		
		newLink.innovation = newInnovation();
		newLink.weight = Math.random() * 4-2;
		
		
		genome.genes.push(newLink);
	}
	
	var nodeMutate = function(genome){
		if(genome.genes.length == 0){
			return;
		}
		
		genome.maxneuron++;
		
		var gene = genome.genes[parseInt(Math.random() * genome.genes.length)];
		if(!gene.enable) return;
		
		gene.enable = false;
		
		var gene1 = copyGene(gene);
		gene1.out = genome.maxneuron;
		gene1.weight = 1.0;
		gene1.innovation = newInnovation();
		gene1.ennable = true;
		genome.genes.push(gene1);
		
		var gene2 = copyGene(gene);
		gene2.into = genome.maxneuron;
		gene2.innovation = newInnovation();
		gene2.enable = true;
		genome.genes.push(gene2);
		console.log('nodeMuteate');
	}
	
	var enableDisableMutate = function(genome,enable) {
		var candidates = [];
		for(var k in genome.genes){
			var gene = genome.genes[k];
			if(gene.enable != enable) {
				candidates.push(gene);
			}
		}
		
		if(candidates.length == 0){
			return;
		}
		
		var gene = candidates[parseInt(Math.random() * candidates.length)];
		gene.enabled = !gene.enabled;
	}
	
	var mutate = function(genome){
		for(var k in genome.mutationRates){
			if(parseInt(Math.random() * 2) == 1) {
				genome.mutationRates[k] *= 0.95;
			} else {
				genome.mutationRates[k] *= 1.05;
			}
		}
		
		if(Math.random() < genome.mutationRates['connections']) {
			pointMutate(genome);
		}
		
		var p = genome.mutationRates['link'];
		while(p>0){
			if(Math.random() < p) {
				linkMutate(genome,false);
			}
			p--;
		}
		
		p = genome.mutationRates['bias'];
		while(p > 0){
			if(Math.random() < p){
				linkMutate(genome,true);
			}
			p--;
		}
		
		p = genome.mutationRates['node'];
		while(p > 0){
			if(Math.random() < p){
			nodeMutate(genome);
			}
			p--;
		}
		
		p = genome.mutationRates['enable'];
		while(p > 0){
			if(Math.random() < p){
				enableDisableMutate(genome,true);
			}
			p--;
		}
		
		p = genome.mutationRates['disable'];
		while(p > 0){
			if(Math.random() < p){
				enableDisableMutate(genome,false);
			}
			p--;
		}
	}
	
	var disjoin = function(genes1, genes2) {
		var i1 = {};
		for(var i = 0; i < genes1.length; i++){
			var gene = genes1[i];
			i1[gene.innovation] = true;
		}
		
		var i2 = {};
		for(var i = 0; i < genes2.lenth; i++){
			var gene = genes2[i];
			i2[gene.innovation] = true;
		}
		
		var disjointGenes = 0;
		for(var i = 1; i < genes1.length; i++){
			var gene = genes1[i];
			if(!i2[gene.innovation]){
				disjointGenes++;
			}
		}
		
		for(var i = 1; i < genes2.length; i++) {
			var gene = genes2[i];
			if(!i1[gene.innovation]){
				disjointGenes++;
			}
		}
		
		var n = Math.max(genes1.length,genes2.length);
		if(n==0) n = 1;
		return disjointGenes / n;
	}
	
	var weights = function(genes1, genes2){
		var i2 = {};
		for(var i = 1; i < genes2.length; i++){
			i2[genes2[i].innovation] = genes2[i];
		}
		
		var sum = 0;
		var coincident = 0;
		
		for(var i = 0; i < genes1.length; i++){
			var gene = genes1[i];
			if(i2[gene.innovation] != null) {
				var gene2 = i2[gene.innovation];
				sum += Math.abs(gene.weight - gene2.weight);
				coincident++;
			}
		}
		if(coincident == 0) coincident=1;
		return sum / coincident;
	}
	
	var sameSpecies = function(genome1, genome2){
		var dd = deltaValues.deltaDisjoint * disjoin(genome1.genes, genome2.genes);
		var dw = deltaValues.deltaWeights * weights(genome1.genes, genome2.genes);
		
		//console.log(dd,dw,deltaValues.deltaThreshold,dd+dw < deltaValues.deltaThreshold);
		return dd+dw < deltaValues.deltaThreshold;
	}
	
	var rankGlobally = function(){
		var global = [];
		for(var k in pool.species){
			var specie = pool.species[k];
			for(var kG in specie.genomes) {
				global.push(specie.genomes[kG]);
			}
		}
		global.sort(function(a,b){return a.fitness - b.fitness});
		
		for(var i = 0; i < global.length; i++){
			global[i].globalRank = i;
		}
	}

	var calculateAverageFitness = function(species){
		var total = 0;
		
		for(var k in species.genomes) {
			var genome = species.genomes[k];
			total += genome.globalRank;
		}
		
		species.averageFitness = total / species.genomes.length;
	}
	
	var totalAverageFitness = function(){
		var total = 0;
		for(var k in pool.species) {
			total += pool.species[k].averageFitness;
		}
		return total;
	}
	
	var cullSpecies = function(cutToOne){
		for(var kS in pool.species){
			var specie = pool.species[kS];
			specie.genomes.sort(function(a,b){return a.fitness - b.fitness});
			
			var remaining = Math.ceil(specie.genomes.length / 2);
			
			if(cutToOne){
				remaining = 1;
			}
			
			while(specie.genomes.length > remaining){
				specie.genomes.shift();
			}
		}
	}
	
	var breedChild = function(specie){
		var child = {};
		if(Math.random() < chances.crossoverChance) {
			g1 = specie.genomes[parseInt(Math.random() * specie.genomes.length)];
			g2 = specie.genomes[parseInt(Math.random() * specie.genomes.length)];
			child = crossover(g1,g2);
		} else {
			g = specie.genomes[parseInt(Math.random() * specie.genomes.length)];
			child = copyGenome(g);
		}
		
		mutate(child);
		return(child);
	}
	
	var removeStaleSpecies = function(){
		var survived = [];
		
		for(var sK in pool.species){
			var specie = pool.species[sK];
			
			specie.genomes.sort(function(a,b){return b.fitness - a.fitness});
			
			if(specie.genomes[0].fitness > specie.topFitness){
				specie.topFitness = specie.genomes[0].fitness;
				specie.staleness = 0;
			} else {
				specie.staleness++;
			}
			
			if(specie.staleness < staleSpecies || specie.topFitness >= pool.maxFitness){
				survived.push(specie);
			}
		}
		
		pool.species = survived;
	}
	
	var removeWeakSpecies = function(){
		var survived = [];
		
		var sum = totalAverageFitness();
		for(var kS in pool.species){
			var specie = pool.species[kS];
			var breed = Math.floor(specie.averageFitness / sum * population);
			if(breed >= 1){
				survived.push(specie);
			}
		}
		
		pool.species = survived;
	}
	
	var addToSpecies = function(child){
		var foundSpecies = false;
		for(var kS in pool.species){
			var specie = pool.species[kS];
			if(sameSpecies(child,specie.genomes[0])) {
				specie.genomes.push(child);
				foundSpecies = true;
				break;
			}
		}
		
		if(!foundSpecies){
			var childSpecie = newSpecie();
			childSpecie.genomes.push(child);
			pool.species.push(childSpecie);
		}
	}
	
	var initializePool = function(){
		pool = newPool();
		
		for(var i = 0; i < population; i++){
			var basic = basicGenome();
			addToSpecies(basic);
		}
	}

	var setMaxValues = function(){
		var topFitness = 0;
		for(var kS in pool.species){
			var specie = pool.species[kS];
			specie.genomes.sort(function(a,b){return a.fitness-b.fitness});
			specie.topFitness = specie.genomes[0].fitness;
			if(topFitness < specie.topFitness){
				topFitness = specie.topFitness;
			}
		}
		pool.maxFitness = topFitness;
	}
		
	initializePool();
	
	return {
		setInputSize: function(size){
			inputSize = size;
		},
		
		setOutputSize: function(size){
			outputSize = size;
		},
		
		getGenomes: function(){
			var genomes = [];
			for(var k in pool.species){
				var specie = pool.species[k];
				for(var kG in specie.genomes){
					genomes.push(specie.genomes[kG]);
				}
			}
			return genomes;
		},
		
		getSummary: function(){
			return "Generation: " + pool.generation
					+ "\nSpecies: " + pool.species.length
					+ "\nInnovation: " + pool.innovation
					+ "\nMaxFitness: " + pool.maxFitness
			;
		},
		
		getGeneration: function(){return pool.generation},
		
		getBetterGenome: function(){
			var genomes = this.getGenomes();
			genomes.sort(function(a,b){return b.fitness - a.fitness;});
			return genomes[0];
		},
		
		newGeneration: function(){
			setMaxValues();
			cullSpecies(false);
			rankGlobally();
			removeStaleSpecies();
			rankGlobally();
			
			for(var kS in pool.species){
				var specie = pool.species[kS];
				calculateAverageFitness(specie);
			}
			
			removeWeakSpecies();
			var sum = totalAverageFitness();
			var children = [];
			for(var kS in pool.species){
				var specie = pool.species[kS];
				var breed = Math.floor(specie.averageFitness / sum * population) -1;
				for(var i = 0; i < breed; i++){
					children.push(breedChild(specie));
				}
			}
			
			cullSpecies(true);
			while(children.length + pool.species.length < population){
				var specie = pool.species[parseInt(Math.random() * pool.species.length)];
				children.push(breedChild(specie));
			}
			
			for(var kC in children){
				var child = children[kC];
				addToSpecies(child);
			}
			
			pool.generation++;
		},
	
		generateNetwork: function(genome){
			var network = {};
			network.neurons = {};
			
			for(var i = 0; i < inputSize; i++){
				network.neurons[i] = newNeuron();
			}
			
			for(var o = 0; o < outputSize; o++){
				network.neurons[maxNodes + o] = newNeuron();
			}
			
			genome.genes.sort(function(a,b){return b.out - a.out});
			
			for(var k in genome.genes){
				var gene = genome.genes[k];
				if(!gene.enabled) continue;
				
				if(network.neurons[gene.out] == null) {
					network.neurons[gene.out] = newNeuron();
				}
				
				var neuron = network.neurons[gene.out];
				neuron.incoming.push(gene);
				
				if(network.neurons[gene.into] == null){
					network.neurons[gene.into] = newNeuron();
				}
			}
			
			genome.network = network;
		},
	
		evaluateGenome: function(genome,inputs) {
			//generateNetwork(genome);
			//console.log(genome);
			output = evaluateNetwork(genome.network, inputs);
			return output;
		},
		
		save: function(){
			localStorage.setItem('pool',JSON.stringify(pool));
		},
		
		load: function(){
			pool = JSON.parse(localStorage.getItem('pool'));
		}
	
	
	}
}
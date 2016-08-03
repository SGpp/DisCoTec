/*
 *  MPI-FT_bitflips.cpp
 *
 *  Created on: 03.01.2016
 *      Author: Johannes Walter
 */
#include <random>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "MPI-FT.h"
#include REAL_MPI_INCLUDE
#define BYTE_SIZE 8

void simft::Sim_FT_Manipulate_bits( int * IntArray, int ArraySize, double p ){
	int S = sizeof(*IntArray); //get the size of one integer

	std::random_device rd; //do we need to seed with time(0)?
	std::mt19937 gen(rd());


	//Calculate amount of bits in the array
	int BitAmount = S * ArraySize * BYTE_SIZE;

	//given the probability p, generate the amount of bits to be flipped randomly using the binomial distribution
	int BitsToFlip = getAmountOfFlips( BitAmount, p );

	std::vector<int> Bits;
	Bits.reserve(BitAmount);

	if(BitsToFlip > 0){

		/*
		 * We want to calculate the exact bit positions in our array to be flipped
		 */

		//fill a vector with numbers 0,1,2,...,BitAmount-1
		for (int i = 0; i < BitAmount; i++) Bits.push_back(i);

		//make sure, the first "BitsToFlip" entries of our vector are randomly selected from the whole vector
		for(int i = 0; i < BitsToFlip; i++){
			int temp;

			//create a random number between the current entry i and the bit amount
			std::uniform_int_distribution<> dis(i, BitAmount-1);
			int randNr = dis(gen);

			//swap the current vector entry with the randomly selected one (also possibly with itself)
			temp = Bits[randNr];
			Bits[randNr] = Bits[i];
			Bits[i] = temp;

			/* Immediately flip the randomly selected bit Bits[i].
			 * Note: alternatively instead of flip it could be possible to always set a defective bit to 1 or always 0 or select it randomly.
			*/
			simft::flipBit(IntArray, Bits[i]);
		}

	}

}

int simft::getAmountOfFlips( int n, double p ){
	std::random_device rd; //do we need to seed with time(0)?
	std::mt19937 gen(rd());
	std::binomial_distribution<> binomial(n, p);
	return binomial(gen);
}

void simft::flipBit (int * BitArray, int BitPos){
	int S = sizeof(*BitArray);
	//calculate the array element, where the bitflip is located
	int ArrayPos = (int) BitPos / (int) (BYTE_SIZE*S);
	//calculate the position of the to-be-flipped bit inside the array element
	int InBitPos = BitPos % (BYTE_SIZE*S);
	//flip the bit. This is achieved by using XOR with 1 on the defective bit
	BitArray[ArrayPos] = BitArray[ArrayPos] ^ (0x1 << ((BYTE_SIZE*S - 1) - InBitPos));;
}

/*
  \file   BIBASPIRfile.h
  \brief  Definitions useful for the BIBASPIR file format.
*/

#ifndef _BIBASPIR_FILE_H_
#define _BIBASPIR_FILE_H_

//#define BIBASPIR_BYTES                320012288
#define BIBASPIR_SAMPLING_RATE         80000000.0
#define BIBASPIR_MSECS_IN_FILE              999
//#define BIBASPIR_WAV_LENGTH                 512
#define BIBASPIR_NBLOCKS                   4882                                //  4.882 number of data blocks
#define BIBASPIR_BLOCK_SIZE            (64*1024)                               // 65.536 size of block (Bytes)
#define BIBASPIR_ID_SIZE                      4                                //      4 size of ID
#define BIBASPIR_PAYLOAD_SIZE      (BIBASPIR_BLOCK_SIZE-BIBASPIR_ID_SIZE)      // 65.532 payload size for 1 block (Bytes)
#define BIBASPIR_PAYLOAD_SAMPLES   (BIBASPIR_PAYLOAD_SIZE/4)                   // 16.383 number of sampling instants in 1 block
#define BIBASPIR_NSAMPLES          (BIBASPIR_PAYLOAD_SAMPLES*BIBASPIR_NBLOCKS) //  79.981.806 sampling instants in one file
#define BIBASPIR_FILESIZE          (BIBASPIR_BLOCK_SIZE*(BIBASPIR_NBLOCKS+1))  // 320.012.288 total length of file
#define BIBASPIR_ANCILLARYPOS      (BIBASPIR_BLOCK_SIZE*BIBASPIR_NBLOCKS)      // 319.946.752 position at which ancillary data starts
#define BIBASPIR_FLAG_POS          (BIBASPIR_ANCILLARYPOS+0xC1)                // 319.946.948 position of the flag that defines the filetype
#define BIBASPIR_FLAG              0x42         // 'B' for BIBASPIR      
#define BIBASPIR_BLOCK_ID_MASK     0x0000FFFF
#define BIBASPIR_ANCI_ID_BIT       30
#define BIBASPIR_ANCI_ID_MASK      0x00000003

#define BIBASPIR_BLOCK_SIZE_INTS   (BIBASPIR_BLOCK_SIZE/4)
#define BIBASPIR_ID_SIZE_INTS      (BIBASPIR_ID_SIZE/4)

#define BIBASPIR_UP_MASK           0x0000FFFF
#define BIBASPIR_L1_MASK           0x00FF00FF
#define BIBASPIR_I_MASK            0xF0F0F0F0

#endif

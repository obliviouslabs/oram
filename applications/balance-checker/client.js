const crypto = require('crypto');
const fetch = require('node-fetch');

// Constants
const IV_SIZE = 12; // AES-GCM IV size

// Utility functions for Base64 encoding/decoding
const base64Encode = (bytes) => Buffer.from(bytes).toString('base64');
const base64Decode = (str) => Uint8Array.from(Buffer.from(str, 'base64'));

// Derive a shared secret
function deriveKey(publicKey) {
    const ecdh = crypto.createECDH('prime256v1');
    ecdh.generateKeys();
    return [ecdh.computeSecret(publicKey), ecdh.getPublicKey()];
}

// Encrypt a message
function encrypt(message, key, iv) {
    const cipher = crypto.createCipheriv('aes-256-gcm', key, iv);
    cipher.setAutoPadding(false);
    const encrypted = Buffer.concat([cipher.update(message), cipher.final()]);
    return {
        encryptedData: encrypted,
        authTag: cipher.getAuthTag()
    };
}

// Decrypt a message
function decrypt(encrypted, key, iv, authTag) {
    const decipher = crypto.createDecipheriv('aes-256-gcm', key, iv);
    decipher.setAutoPadding(true);
    decipher.setAuthTag(authTag);
    const decrypted = Buffer.concat([decipher.update(encrypted), decipher.final()]);
    return decrypted;
}

// example balance checker client
(async function () {
    const host = "http://localhost:8080";
    const addr = "0xa7C0D36c4698981FAb42a7d8c783674c6Fe2592d";
    try {
        const res = await fetch(`${host}/public_key`);
        if (!res.ok) throw new Error(`HTTP error! status: ${res.status}`);

        const serverPubKeyBase64 = await res.text();
        console.log("Response:", serverPubKeyBase64);
        const serverPubKeyBytes = base64Decode(serverPubKeyBase64);

        console.log("Server public key received");

        // // print the first 32 bytes in hex
        const [sharedSecret, clientPublicKey] = deriveKey(Buffer.concat([Buffer.from([0x04]), serverPubKeyBytes]));

        const clientPubKeyBytes = new Uint8Array(clientPublicKey);
        console.log("Shared secret computed");

        nonce = crypto.randomBytes(4)
        const body = await makeBalanceQueryBody('USDT', addr, nonce, clientPubKeyBytes, sharedSecret);

        const postRes = await fetch(`${host}/secure`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/octet-stream',
            },
            body: body,
        });


        const encryptedResponse = await postRes.text();

        const { balance, lastBlock, successFlag, queryTime } = await decodeResponse(encryptedResponse, nonce, sharedSecret);
        console.log("Balance:", balance);
        console.log("Last Block:", lastBlock);
        console.log("Success Flag:", successFlag);
        console.log("Query Time:", queryTime);

    } catch (error) {
        console.error("Error:", error);
    }

    async function makeBalanceQueryBody(coinType, addr, nonce, clientPubKey, sharedSecret) {
        // hex to bytes
        if (addr.length !== 42 || !addr.startsWith("0x")) {
            throw new Error("Invalid address format");
        }
        const addrBytes = Buffer.from(addr.slice(2), 'hex');
        READ_BALANCE = Buffer.alloc(4);
        USDT = Buffer.alloc(4);
        const queryBuf = Buffer.concat([READ_BALANCE, USDT, addrBytes, nonce]);

        const iv = crypto.randomBytes(IV_SIZE);
        const { encryptedData, authTag } = encrypt(queryBuf, sharedSecret, iv);

        console.log("Encrypted Data:", Buffer.from(encryptedData).toString('hex'));
        console.log("Auth Tag:", Buffer.from(authTag).toString('hex'));
        console.log("IV:", Buffer.from(iv).toString('hex'));
        return base64Encode(Buffer.concat([encryptedData, authTag, clientPubKey.slice(1), iv]));
    }

    async function decodeResponse(responseBase64, nonce, sharedSecret) {
        const responseBytes = base64Decode(responseBase64);

        const balanceLen = 32
        const lastBlockLen = 8
        const nonceLen = 4
        const successFlagLen = 4
        const queryTimeLen = 8
        const tagLen = 16
        const ivLen = 12
        const encryptedDataLen = balanceLen + lastBlockLen + nonceLen + successFlagLen + queryTimeLen
        if (responseBytes.length !== (encryptedDataLen + tagLen + ivLen - 1) / 8 * 8 + 1) {
            throw new Error("Invalid encrypted response size, expected 68 bytes, got " + responseBytes.length + " bytes");
        }
        offset = 0
        const encryptedData = responseBytes.slice(offset, offset + encryptedDataLen); // Assuming first IV_SIZE bytes are IV
        offset += encryptedDataLen
        const authTag = responseBytes.slice(offset, offset + tagLen); // Assuming authTag 
        offset += tagLen
        const iv = responseBytes.slice(offset, offset + ivLen); // Assuming last 12 bytes are IV
        console.log("Shared Secret:", Buffer.from(sharedSecret).toString('hex'));
        console.log("Encrypted Data:", Buffer.from(encryptedData).toString('hex'));
        console.log("Auth Tag:", Buffer.from(authTag).toString('hex'));
        console.log("IV:", Buffer.from(iv).toString('hex'));
        const decryptedResponse = decrypt(encryptedData, sharedSecret, iv, authTag);
        offset = 0
        const balanceBytes = decryptedResponse.slice(0, balanceLen);
        offset += balanceLen
        const lastBlockBytes = decryptedResponse.slice(offset, offset + lastBlockLen);
        offset += lastBlockLen
        const responseNonce = decryptedResponse.slice(offset, offset + nonceLen);
        offset += nonceLen
        const successFlagBytes = decryptedResponse.slice(offset, offset + successFlagLen);
        offset += successFlagLen
        const queryTimeBytes = decryptedResponse.slice(offset, offset + queryTimeLen);
        for (let i = 0; i < nonce.length; i++) {
            if (nonce[i] !== responseNonce[i]) {
                throw new Error("Invalid nonce in response");
            }
        }
        const balance = BigInt("0x" + balanceBytes.toString('hex'));
        const lastBlock = BigInt("0x" + lastBlockBytes.toString('hex'));
        const successFlag = BigInt("0x" + successFlagBytes.toString('hex'));
        const queryTime = BigInt("0x" + queryTimeBytes.toString('hex'));
        return { balance, lastBlock, successFlag, queryTime };

    }
})();
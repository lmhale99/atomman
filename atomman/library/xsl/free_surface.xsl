<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" 
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
  xmlns="http://www.w3.org/TR/xhtml1/strict">
  <xsl:output method="html" encoding="utf-8" indent="yes" />
  
  <xsl:template match="free-surface">
    <div>

      <style>
        .aslist {list-style-type: circle; list-style-position: inside; margin: 10px;}
      </style>

      <h1>Free surface parameter set</h1>
      <ul>
        <li><b><xsl:text>ID: </xsl:text></b><xsl:value-of select="id"/></li>
        <li><b><xsl:text>UUID4: </xsl:text></b><xsl:value-of select="key"/></li>
        <li><b><xsl:text>Family: </xsl:text></b>
          <xsl:choose>
            <xsl:when test="system-family-URL">
              <a href="{system-family-URL}"><xsl:value-of select="system-family"/></a>
            </xsl:when>
            <xsl:otherwise>
              <xsl:value-of select="system-family"/>
            </xsl:otherwise>
          </xsl:choose>
        </li>
        <li>
          <b><xsl:text>atomman FreeSurface defect parameters:</xsl:text></b>
          <ul class="aslist">
            <li><xsl:text>hkl: "</xsl:text><xsl:value-of select="calculation-parameter/hkl"/><xsl:text>"</xsl:text></li>
            <xsl:if test="calculation-parameter/shiftindex">
              <li><xsl:text>shiftindex: "</xsl:text><xsl:value-of select="calculation-parameter/shiftindex"/><xsl:text>"</xsl:text></li>
            </xsl:if>
            <li><xsl:text>cutboxvector: "</xsl:text><xsl:value-of select="calculation-parameter/cutboxvector"/><xsl:text>"</xsl:text></li>
          </ul>
        </li>
      </ul>
    </div>

  </xsl:template>
</xsl:stylesheet>